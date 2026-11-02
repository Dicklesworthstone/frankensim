//! Exactly incompressible plane stress from the original matrix/fibre energy.
//!
//! With `a = [F00,F01,F10,F11]` and `d = det(a) > 0`, set `F22 = 1/d`
//! and relax `s = [F02,F12]` in `F = [a0,a1,s0; a2,a3,s1; 0,0,1/d]`.
//! Thus `J=1` is kinematic, independent of every bulk parameter. The energy
//! is [`Hyperelastic::isochoric_energy`], including its original NH remainder
//! and full-invariant HGO fibres. No bulk penalty is evaluated.
//!
//! The pressure multiplier satisfies `P = dW_remainder/dF - p cof(F)`.
//! Eliminating `P22` gives `p = (dW_remainder/dF22)/d`; the remaining two
//! transverse tractions are the shear stationarity equations. Different
//! off-constraint energy extensions shift this multiplier, so it is not the
//! negative mean Cauchy stress. For the native NH remainder, `p=mu*(F22²-1)`.
//!
//! The tangent differentiates the six-variable constrained energy, including
//! the first AND second derivatives of `1/d`, then Schur-condenses the two
//! relaxed shears. The incompressible thickness elimination is described by
//! Ahmadi et al. (2025), doi:10.1007/s00707-025-04310-z; the transverse-shear
//! extension here follows directly from stationarity for tilted fibres.
//! This explicit constitutive choice is a local stable membrane branch, not
//! a 3-D mixed formulation, shell, wrinkling or global-uniqueness certificate.

use super::{Hyperelastic, plane_stress::PlaneStressSettings};
use crate::MaterialError;
use fs_ad::{Dual, Dual64, Real};

const IN_PLANE: [usize; 4] = [0, 1, 3, 4];
const SHEAR: [usize; 2] = [2, 5];

/// Converged exactly incompressible plane-stress energy and stress.
#[derive(Clone, Copy, Debug)]
pub struct IncompressiblePlaneStressState {
    /// Full row-major deformation gradient, with `F22=1/det(F_in_plane)`.
    pub deformation_gradient: [f64; 9],
    /// Original remainder energy per reference volume (Pa = J/m³).
    pub energy_density_pa: f64,
    /// In-plane first Piola stress `[P00,P01,P10,P11]`, in Pa.
    pub piola: [f64; 4],
    /// Full physical first Piola stress, including the pressure multiplier.
    pub full_piola: [f64; 9],
    /// Compression-positive multiplier in `P = dW_remainder/dF - p cof(F)`.
    /// Its value depends on the native remainder's pressure convention.
    pub pressure_pa: f64,
    /// Transverse traction norm divided by the fixed initial shear scale.
    pub relative_transverse_residual: f64,
    /// Number of accepted transverse-shear Newton updates.
    pub iterations: usize,
}

/// Exactly incompressible response with its consistently reduced tangent.
#[derive(Clone, Copy, Debug)]
pub struct IncompressiblePlaneStressResponse {
    /// Full row-major deformation gradient, with `F22=1/det(F_in_plane)`.
    pub deformation_gradient: [f64; 9],
    /// Original remainder energy per reference volume (Pa = J/m³).
    pub energy_density_pa: f64,
    /// In-plane first Piola stress `[P00,P01,P10,P11]`, in Pa.
    pub piola: [f64; 4],
    /// Full physical first Piola stress, including the pressure multiplier.
    pub full_piola: [f64; 9],
    /// Compression-positive multiplier for the native energy remainder.
    pub pressure_pa: f64,
    /// `dP_in_plane/dF_in_plane`, including thickness and shear relaxation.
    pub tangent: [[f64; 4]; 4],
    /// Transverse traction norm divided by the fixed initial shear scale.
    pub relative_transverse_residual: f64,
    /// Number of accepted transverse-shear Newton updates.
    pub iterations: usize,
}

fn state_error(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("incompressible hyperelastic plane stress: {what}"),
    }
}

// Six independent coordinates. The determinant constraint is differentiated
// through this map, not imposed after differentiating a compressible tangent.
fn constrained_gradient<T: Real>(x: &[T; 6]) -> [T; 9] {
    let thickness = T::one() / (x[0] * x[3] - x[1] * x[2]);
    [
        x[0],
        x[1],
        x[4],
        x[2],
        x[3],
        x[5],
        T::zero(),
        T::zero(),
        thickness,
    ]
}

// Scalar-normalized Cholesky of the two eliminated shear coordinates.
// An unstable block refuses; no regularization or stiffness floor is used.
struct PositiveShear {
    lower: [[f64; 2]; 2],
    scale: f64,
}

impl PositiveShear {
    fn new(a: [[f64; 2]; 2]) -> Result<Self, MaterialError> {
        let scale = a.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
        if !(scale.is_finite() && scale > 0.0)
            || a.iter().flatten().any(|v| !v.is_finite())
            || (a[0][1] / scale - a[1][0] / scale).abs() > 128.0 * f64::EPSILON
        {
            return Err(state_error(
                "nonfinite, zero or nonsymmetric shear stiffness",
            ));
        }
        let a00 = a[0][0] / scale;
        if a00 <= 0.0 {
            return Err(state_error("transverse shear branch is not stable"));
        }
        let l00 = fs_math::det::sqrt(a00);
        let l10 = (a[1][0] / scale) / l00;
        let remainder = a[1][1] / scale - l10 * l10;
        if !(remainder.is_finite() && remainder > 0.0) {
            return Err(state_error("transverse shear branch is not stable"));
        }
        Ok(Self {
            lower: [[l00, 0.0], [l10, fs_math::det::sqrt(remainder)]],
            scale,
        })
    }

    fn solve(&self, rhs: [f64; 2]) -> Result<[f64; 2], MaterialError> {
        let l = &self.lower;
        let y0 = (rhs[0] / self.scale) / l[0][0];
        let y1 = (rhs[1] / self.scale - l[1][0] * y0) / l[1][1];
        let x1 = y1 / l[1][1];
        let x0 = (y0 - l[1][0] * x1) / l[0][0];
        if !x0.is_finite() || !x1.is_finite() {
            return Err(state_error("unrepresentable transverse shear solve"));
        }
        Ok([x0, x1])
    }

    fn traction_scale(&self) -> Result<f64, MaterialError> {
        let compliance = self.solve([self.scale, 0.0])?[0] + self.solve([0.0, self.scale])?[1];
        let scale = self.scale / compliance;
        if !(scale.is_finite() && scale > 0.0) {
            return Err(state_error("unrepresentable shear compliance scale"));
        }
        Ok(scale)
    }
}

impl Hyperelastic {
    fn incompressible_shear_data(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 2], PositiveShear), MaterialError> {
        self.check_state(f)?;
        let mut stress = [0.0; 2];
        let mut hessian = [[0.0; 2]; 2];
        for col in 0..2 {
            let mut vars = f.map(|v| Dual::<Dual64<2>, 1>::constant(Dual64::constant(v)));
            for (row, &index) in SHEAR.iter().enumerate() {
                vars[index] = Dual {
                    re: Dual64::variable(f[index], row),
                    eps: [Dual64::constant(if row == col { 1.0 } else { 0.0 })],
                };
            }
            let w = self.isochoric_energy(&vars);
            if !w.re.re.is_finite() {
                return Err(state_error("nonfinite remainder energy"));
            }
            stress = w.re.eps;
            for row in 0..2 {
                hessian[row][col] = w.eps[0].eps[row];
            }
        }
        if stress.iter().any(|v| !v.is_finite()) {
            return Err(state_error("nonfinite transverse shear stress"));
        }
        Ok((stress, PositiveShear::new(hessian)?))
    }

    /// Choose exact incompressibility, impose `F22=1/det(F_in_plane)`, and
    /// relax the two transverse shears without evaluating any bulk penalty.
    ///
    /// The initial shears are zero. The original principal-stretch domain and
    /// positive shear stiffness are required at every accepted iterate.
    /// Newton updates and residual-decreasing backtracks use the supplied
    /// budgets. The relative tolerance uses `1/trace(Hss^-1)` times
    /// `max(1, ||F_in_plane||)` at the initial state; changing stress units
    /// cannot change the selected equilibrium. Callers poll cancellation
    /// between quadrature points. The finite bulk parameter remains stored in
    /// the law but is deliberately irrelevant to this explicit material mode.
    pub fn incompressible_plane_stress_state(
        &self,
        in_plane: &[f64; 4],
        settings: PlaneStressSettings,
    ) -> Result<IncompressiblePlaneStressState, MaterialError> {
        if !(settings.relative_tolerance.is_finite()
            && settings.relative_tolerance > 0.0
            && settings.relative_tolerance < 1.0)
            || settings.max_iterations == 0
            || settings.max_backtracks == 0
        {
            return Err(MaterialError::Parameters {
                what: "incompressible plane stress needs 0 < relative_tolerance < 1 and positive iteration/backtracking budgets".into(),
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
        let mut f =
            constrained_gradient(&[in_plane[0], in_plane[1], in_plane[2], in_plane[3], 0.0, 0.0]);
        let (mut stress, mut block) = self.incompressible_shear_data(&f)?;
        let gradient_norm = in_plane[0]
            .hypot(in_plane[1])
            .hypot(in_plane[2])
            .hypot(in_plane[3]);
        let scale = block.traction_scale()? * gradient_norm.max(1.0);
        if !(scale.is_finite() && scale > 0.0) {
            return Err(state_error("unrepresentable traction normalization"));
        }
        for iteration in 0..=settings.max_iterations {
            let residual = (stress[0] / scale).hypot(stress[1] / scale);
            if residual <= settings.relative_tolerance {
                let (energy, remainder_piola) =
                    fs_ad::gradient::<9>(f, |x| self.isochoric_energy(&x));
                let pressure = remainder_piola[8] / determinant;
                let cofactor = [
                    f[4] * f[8],
                    -f[3] * f[8],
                    0.0,
                    -f[1] * f[8],
                    f[0] * f[8],
                    0.0,
                    f[1] * f[5] - f[2] * f[4],
                    f[2] * f[3] - f[0] * f[5],
                    determinant,
                ];
                let full_piola =
                    core::array::from_fn(|i| remainder_piola[i] - pressure * cofactor[i]);
                let final_residual = (full_piola[2] / scale)
                    .hypot(full_piola[5] / scale)
                    .hypot(full_piola[8] / scale);
                if !energy.is_finite()
                    || !pressure.is_finite()
                    || full_piola.iter().any(|v| !v.is_finite())
                    || !(final_residual <= settings.relative_tolerance)
                {
                    return Err(state_error(
                        "unrepresentable or unresolved physical plane stress",
                    ));
                }
                return Ok(IncompressiblePlaneStressState {
                    deformation_gradient: f,
                    energy_density_pa: energy,
                    piola: IN_PLANE.map(|i| full_piola[i]),
                    full_piola,
                    pressure_pa: pressure,
                    relative_transverse_residual: final_residual,
                    iterations: iteration,
                });
            }
            if iteration == settings.max_iterations || !residual.is_finite() {
                break;
            }
            let direction = block.solve(stress.map(|v| -v))?;
            let mut alpha = 1.0;
            let mut accepted = None;
            for _ in 0..settings.max_backtracks {
                let mut trial = f;
                for i in 0..2 {
                    trial[SHEAR[i]] += alpha * direction[i];
                }
                match self.incompressible_shear_data(&trial) {
                    Ok((next_stress, next_block)) => {
                        let next_residual = (next_stress[0] / scale).hypot(next_stress[1] / scale);
                        if next_residual <= (1.0 - 1e-4 * alpha) * residual {
                            accepted = Some((trial, next_stress, next_block));
                            break;
                        }
                    }
                    Err(MaterialError::State { .. }) => (),
                    Err(error) => return Err(error),
                }
                alpha *= 0.5;
            }
            let Some((next_f, next_stress, next_block)) = accepted else {
                return Err(state_error(
                    "transverse shear line search exhausted its budget or material domain",
                ));
            };
            f = next_f;
            stress = next_stress;
            block = next_block;
        }
        Err(state_error(
            "transverse shear Newton solve exhausted its iteration budget",
        ))
    }

    /// Exact-incompressibility state plus the constrained-energy tangent,
    /// including consistent elimination of both transverse shear coordinates.
    /// A six-direction AD Hessian is formed only after local equilibrium.
    pub fn incompressible_plane_stress(
        &self,
        in_plane: &[f64; 4],
        settings: PlaneStressSettings,
    ) -> Result<IncompressiblePlaneStressResponse, MaterialError> {
        let state = self.incompressible_plane_stress_state(in_plane, settings)?;
        let f = state.deformation_gradient;
        let x = [
            in_plane[0],
            in_plane[1],
            in_plane[2],
            in_plane[3],
            f[2],
            f[5],
        ];
        let mut hessian = [[0.0; 6]; 6];
        for col in 0..6 {
            let vars = core::array::from_fn(|row| Dual {
                re: Dual64::<6>::variable(x[row], row),
                eps: [Dual64::constant(if row == col { 1.0 } else { 0.0 })],
            });
            let w = self.isochoric_energy(&constrained_gradient(&vars));
            for row in 0..6 {
                hessian[row][col] = w.eps[0].eps[row];
            }
        }
        if hessian.iter().flatten().any(|v| !v.is_finite()) {
            return Err(state_error("nonfinite constrained-energy Hessian"));
        }
        let block = PositiveShear::new([
            [hessian[4][4], hessian[4][5]],
            [hessian[5][4], hessian[5][5]],
        ])?;
        let mut tangent = [[0.0; 4]; 4];
        for col in 0..4 {
            let extension = block.solve([hessian[4][col], hessian[5][col]])?;
            for row in 0..4 {
                tangent[row][col] = hessian[row][col]
                    - hessian[row][4] * extension[0]
                    - hessian[row][5] * extension[1];
            }
        }
        if tangent.iter().flatten().any(|v| !v.is_finite()) {
            return Err(state_error("unrepresentable condensed tangent"));
        }
        Ok(IncompressiblePlaneStressResponse {
            deformation_gradient: f,
            energy_density_pa: state.energy_density_pa,
            piola: state.piola,
            full_piola: state.full_piola,
            pressure_pa: state.pressure_pa,
            tangent,
            relative_transverse_residual: state.relative_transverse_residual,
            iterations: state.iterations,
        })
    }
}
