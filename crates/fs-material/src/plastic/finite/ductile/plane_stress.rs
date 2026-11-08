//! Traction-free finite sheets with progressive ductile softening.
//!
//! The prescribed in-plane gradient is `[F00,F01,F10,F11]`. All three entries
//! of the transverse column `[F02,F12,F22]` are solved from `P[:,2]=0`, with
//! positive thickness. The complete original ductile return, accepted history,
//! temperature and physical interval stay fixed throughout the local solve.
//!
//! Damage makes the discrete tangent generally nonsymmetric. A pivoted general
//! transverse solve therefore supplies both Newton updates and the exact Schur
//! complement `Aaa - Aat Att^-1 Ata`; no symmetrization, stiffness floor or
//! scalar incremental potential is introduced. This selects a regular local
//! equilibrium branch. It does not certify its stability or global uniqueness.
//! Shell bending, necking, tool contact and crack topology remain separate.
//! Exact finite-support failure removes the transverse shear constraints.
//! Terminal candidates refuse before dividing by integrity or inverting the
//! transverse block; no arbitrary thickness, stiffness floor or deletion is
//! substituted for the lost three-dimensional material constraint.

use super::{FiniteDuctileError, FiniteDuctileJ2, FiniteDuctileState, FiniteDuctileStep};
use crate::ductile::evolution::{DifferentiableDuctileInitiationLaw, DuctileJ2Conditions};
use crate::hyper::plane_stress::PlaneStressSettings;
use crate::plastic::finite::FiniteJ2Error;

const IN_PLANE: [usize; 4] = [0, 1, 3, 4];
const TRANSVERSE: [usize; 3] = [2, 5, 8];

/// Disposable plane-stress proposal, including the original full work account.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteDuctilePlaneStressStep<'a, L> {
    /// Full row-major gradient, including solved transverse shear and thickness.
    pub deformation_gradient: [f64; 9],
    /// Physical in-plane first Piola stress `[P00,P01,P10,P11]`, Pa.
    pub piola_pa: [f64; 4],
    /// Exact branch derivative after eliminating transverse equilibrium.
    /// This can be nonsymmetric and indefinite during progressive softening.
    pub tangent_pa: [[f64; 4]; 4],
    /// Original full candidate, stress, history and reference-volume energy.
    /// Storage, plastic loss, released damage storage and signed endpoint-work
    /// defect retain their native meanings. No extra loss is created here.
    pub step: FiniteDuctileStep<'a, L>,
    /// Norm of transverse traction divided by surviving integrity and a fixed
    /// initial compliance scale. Vanishing integrity cannot hide a bad root.
    pub relative_transverse_residual: f64,
    /// Accepted local updates, excluding initialization and backtracks.
    pub iterations: usize,
}

impl<L> FiniteDuctilePlaneStressStep<'_, L> {
    /// Current normal thickness per unit original thickness.
    pub fn thickness_stretch(&self) -> f64 {
        self.deformation_gradient[8]
    }

    /// Warm start for a later endpoint; this does not replace material history.
    pub fn transverse_column(&self) -> [f64; 3] {
        TRANSVERSE.map(|i| self.deformation_gradient[i])
    }
}

fn numerical(what: &'static str) -> FiniteDuctileError {
    FiniteJ2Error::Numerical(what).into()
}

// Scalar-normalized partial-pivot LU on the actual nonsymmetric 3x3 block.
// Fixed-size storage avoids an allocation for every material point and retry.
struct TransverseBlock {
    factors: [[f64; 3]; 3],
    permutation: [usize; 3],
    scale: f64,
}

impl TransverseBlock {
    fn new(tangent: &[[f64; 9]; 9]) -> Result<Self, FiniteDuctileError> {
        let a = TRANSVERSE.map(|i| TRANSVERSE.map(|j| tangent[i][j]));
        let scale = a.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
        if !scale.is_finite() || scale <= 0.0 || a.iter().flatten().any(|v| !v.is_finite()) {
            return Err(numerical(
                "ductile plane-stress transverse stiffness is nonfinite or zero",
            ));
        }
        let mut factors = a.map(|row| row.map(|v| v / scale));
        let mut permutation = [0, 1, 2];
        for column in 0..3 {
            let mut pivot = column;
            for row in column + 1..3 {
                // Strict comparison retains the first row on an exact tie.
                if factors[row][column].abs() > factors[pivot][column].abs() {
                    pivot = row;
                }
            }
            if factors[pivot][column] == 0.0 {
                return Err(numerical(
                    "ductile plane-stress transverse tangent is singular",
                ));
            }
            factors.swap(column, pivot);
            permutation.swap(column, pivot);
            for row in column + 1..3 {
                factors[row][column] /= factors[column][column];
                for j in column + 1..3 {
                    factors[row][j] =
                        (-factors[row][column]).mul_add(factors[column][j], factors[row][j]);
                }
            }
        }
        if factors.iter().flatten().any(|v| !v.is_finite()) {
            return Err(numerical(
                "ductile plane-stress factorization is not representable",
            ));
        }
        Ok(Self {
            factors,
            permutation,
            scale,
        })
    }

    fn solve(&self, rhs: [f64; 3]) -> Result<[f64; 3], FiniteDuctileError> {
        let mut x = self.permutation.map(|i| rhs[i] / self.scale);
        for row in 0..3 {
            for j in 0..row {
                x[row] = (-self.factors[row][j]).mul_add(x[j], x[row]);
            }
        }
        for row in (0..3).rev() {
            for j in row + 1..3 {
                x[row] = (-self.factors[row][j]).mul_add(x[j], x[row]);
            }
            x[row] /= self.factors[row][row];
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(numerical(
                "ductile plane-stress transverse solve is not representable",
            ));
        }
        Ok(x)
    }

    fn traction_scale(&self) -> Result<f64, FiniteDuctileError> {
        let mut compliance_norm = 0.0_f64;
        for j in 0..3 {
            let mut rhs = [0.0; 3];
            rhs[j] = self.scale;
            for value in self.solve(rhs)? {
                compliance_norm = compliance_norm.hypot(value);
            }
        }
        // 1/||Att^-1||_F <= the smallest singular value. A large bulk modulus
        // therefore cannot hide a transverse shear residual, even without
        // symmetry or a positive-definite algorithmic tangent.
        let scale = self.scale / compliance_norm;
        if !scale.is_finite() || scale <= 0.0 {
            return Err(numerical(
                "ductile plane-stress compliance scale is not representable",
            ));
        }
        Ok(scale)
    }
}

fn residual<L>(step: &FiniteDuctileStep<'_, L>, scale: f64) -> f64 {
    let integrity = step.state.integrity();
    TRANSVERSE.iter().fold(0.0_f64, |norm, &i| {
        norm.hypot((step.first_piola_pa[i] / integrity) / scale)
    })
}

fn condense<'a, L>(
    f: [f64; 9],
    step: FiniteDuctileStep<'a, L>,
    block: &TransverseBlock,
    relative_transverse_residual: f64,
    iterations: usize,
) -> Result<FiniteDuctilePlaneStressStep<'a, L>, FiniteDuctileError> {
    let a = &step.tangent_pa;
    let mut tangent_pa = [[0.0; 4]; 4];
    for (j, &col) in IN_PLANE.iter().enumerate() {
        let extension = block.solve(TRANSVERSE.map(|k| a[k][col]))?;
        for (i, &row) in IN_PLANE.iter().enumerate() {
            tangent_pa[i][j] = a[row][col]
                - (0..3)
                    .map(|k| a[row][TRANSVERSE[k]] * extension[k])
                    .sum::<f64>();
        }
    }
    if tangent_pa.iter().flatten().any(|v| !v.is_finite()) {
        return Err(numerical(
            "ductile plane-stress condensed tangent is not representable",
        ));
    }
    Ok(FiniteDuctilePlaneStressStep {
        deformation_gradient: f,
        piola_pa: IN_PLANE.map(|i| step.first_piola_pa[i]),
        tangent_pa,
        step,
        relative_transverse_residual,
        iterations,
    })
}

impl<L: DifferentiableDuctileInitiationLaw + Clone> FiniteDuctileJ2<L> {
    /// Relax all transverse tractions using the full progressive ductile law.
    ///
    /// Every trial borrows the same old state and physical interval, including
    /// selected Perzyna relaxation. Only the final returned candidate may be
    /// committed. Refusal leaves caller history unchanged. Iterations use
    /// fixed-size storage; spatial callers poll cancellation between points.
    ///
    /// The relative tolerance uses `||P[:,2]/g||` and the initial inverse-
    /// compliance scale divided by initial integrity, times `max(1,||F_in||)`.
    /// This keeps traction equilibrium meaningful as the sheet loses strength.
    /// Bounded backtracking preserves positive thickness and original material
    /// limits. A unit-volume seed is only a fallback, never a volume constraint.
    pub fn plane_stress_trial<'a>(
        &'a self,
        in_plane: &[f64; 4],
        old: &FiniteDuctileState<'a, L>,
        transverse_seed: [f64; 3],
        conditions: DuctileJ2Conditions,
        settings: PlaneStressSettings,
    ) -> Result<FiniteDuctilePlaneStressStep<'a, L>, FiniteDuctileError> {
        self.admit(old, conditions)?;
        if old.is_fully_failed() {
            return Err(FiniteDuctileError::PlaneStressTerminalFailure);
        }
        if !(settings.relative_tolerance.is_finite()
            && settings.relative_tolerance > 0.0
            && settings.relative_tolerance < 1.0)
            || settings.max_iterations == 0
            || settings.max_backtracks == 0
        {
            return Err(FiniteDuctileError::Invalid(
                "ductile plane stress requires 0 < relative tolerance < 1 and positive local budgets",
            ));
        }
        let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
        if in_plane
            .iter()
            .chain(&transverse_seed)
            .any(|v| !v.is_finite())
            || !determinant.is_finite()
            || determinant <= 0.0
            || transverse_seed[2] <= 0.0
        {
            return Err(FiniteDuctileError::Invalid(
                "ductile plane stress requires finite orientation-preserving in-plane deformation and positive thickness",
            ));
        }
        let mut f = [
            in_plane[0],
            in_plane[1],
            transverse_seed[0],
            in_plane[2],
            in_plane[3],
            transverse_seed[1],
            0.0,
            0.0,
            transverse_seed[2],
        ];
        let evaluate = |gradient: &[f64; 9]| {
            let step = self.trial(gradient, old, conditions)?;
            if step.state.is_fully_failed() {
                return Err(FiniteDuctileError::PlaneStressTerminalFailure);
            }
            let block = TransverseBlock::new(&step.tangent_pa)?;
            Ok::<_, FiniteDuctileError>((step, block))
        };
        let target_thickness = 1.0 / determinant;
        let (mut step, mut block) = match evaluate(&f) {
            Ok(candidate) => candidate,
            Err(initial_error) => {
                let mut seed = None;
                let mut fraction = 1.0;
                for _ in 0..settings.max_backtracks {
                    let mut trial = f;
                    trial[8] = (1.0 - fraction) * transverse_seed[2] + fraction * target_thickness;
                    if trial[8].is_finite()
                        && trial[8] > 0.0
                        && let Ok((value, factor)) = evaluate(&trial)
                    {
                        seed = Some((trial, value, factor));
                        break;
                    }
                    fraction *= 0.5;
                }
                let Some((gradient, value, factor)) = seed else {
                    return Err(initial_error);
                };
                f = gradient;
                (value, factor)
            }
        };
        let gradient_norm = in_plane.iter().fold(0.0_f64, |n, v| n.hypot(*v));
        let scale = (block.traction_scale()? / step.state.integrity()) * gradient_norm.max(1.0);
        if !scale.is_finite() || scale <= 0.0 {
            return Err(numerical(
                "ductile plane-stress traction normalization is not representable",
            ));
        }
        for iteration in 0..=settings.max_iterations {
            let norm = residual(&step, scale);
            if norm <= settings.relative_tolerance {
                return condense(f, step, &block, norm, iteration);
            }
            if iteration == settings.max_iterations || !norm.is_finite() {
                break;
            }
            let direction = block.solve(TRANSVERSE.map(|i| -step.first_piola_pa[i]))?;
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..settings.max_backtracks {
                let mut trial = f;
                for i in 0..3 {
                    trial[TRANSVERSE[i]] += fraction * direction[i];
                }
                if trial[8] > 0.0
                    && let Ok((value, factor)) = evaluate(&trial)
                    && residual(&value, scale) <= (1.0 - 1e-4 * fraction) * norm
                {
                    accepted = Some((trial, value, factor));
                    break;
                }
                fraction *= 0.5;
            }
            if accepted.is_none() {
                // Strong softening can make the physical-traction Newton
                // direction ascend the integrity-normalized residual far
                // from the root. Try the same bounded volume-neutral seed
                // before refusing; this never advances the old history.
                let mut trial = f;
                trial[8] = target_thickness;
                if trial[8].is_finite()
                    && trial[8] > 0.0
                    && let Ok((value, factor)) = evaluate(&trial)
                    && residual(&value, scale) < norm
                {
                    accepted = Some((trial, value, factor));
                }
            }
            let Some((gradient, value, factor)) = accepted else {
                return Err(numerical(
                    "ductile plane-stress line search exhausted its budget or material domain",
                ));
            };
            f = gradient;
            step = value;
            block = factor;
        }
        Err(numerical(
            "ductile plane-stress Newton solve exhausted its iteration budget",
        ))
    }
}
