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

use super::{Hyperelastic, cards::FiniteStrainPotential};
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
    /// Maximum candidate seeds if unit thickness is unstable, and maximum
    /// trial step lengths per Newton update (1, 1/2, 1/4, ...).
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
    /// Reduced potential per original reference volume (Pa = J/m³).
    /// This is stored energy for elasticity, but an incremental potential for
    /// a frozen inelastic owner; the latter must not be counted as heat.
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
    /// Reduced potential per original reference volume (Pa = J/m³).
    /// Frozen inelastic potentials are not physical stored energy or heat.
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
pub(crate) struct PositiveBlock {
    lower: [[f64; 3]; 3],
    scale: f64,
}
impl PositiveBlock {
    pub(crate) fn new(a: [[f64; 3]; 3]) -> Result<Self, MaterialError> {
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

    pub(crate) fn solve(&self, rhs: [f64; 3]) -> Result<[f64; 3], MaterialError> {
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

    pub(crate) fn traction_scale(&self) -> Result<f64, MaterialError> {
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
}

/// Relax the transverse tractions of the supplied potential without extracting
/// its underlying elastic law. Prescribed reference maps and frozen material
/// histories therefore remain part of every trial and accepted response.
///
/// The original material owns domain admission and the consistent tangent.
/// A frozen inelastic owner must keep its accepted history and physical time
/// interval fixed throughout the solve, and commit its history only afterward.
/// Like the native solve, this selects a local stable branch from unit thickness
/// with bounded volume-preserving fallback seeds; failure to find an admitted
/// seed is a refusal, not a global nonexistence claim.
pub fn plane_stress_state<L: FiniteStrainPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    settings: PlaneStressSettings,
) -> Result<PlaneStressState, MaterialError> {
    generic_state(material, in_plane, None, settings)
}

/// Relax from an explicitly supplied transverse column `[F02,F12,F22]`.
///
/// Prescribed natural configurations can place unit total thickness outside
/// the elastic domain even at zero elastic strain. Their owner should supply
/// an admitted seed, for example `[0,0,det(G)/det(a)]` when zero transverse
/// shear is admitted. This is an initial guess, not a volume constraint: all
/// three entries remain free during the finite-compliance solve.
/// The seed must be finite with positive thickness. If its material state or
/// transverse branch is inadmissible, initialization tries the original unit
/// thickness path before Newton begins. This also retains admissible highly
/// compressible states whose volume-preserving guess leaves the elastic domain.
/// At most one supplied seed precedes the original bounded initialization;
/// no deformation, stiffness or material domain is clipped.
pub fn plane_stress_state_with_seed<L: FiniteStrainPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    initial_transverse: [f64; 3],
    settings: PlaneStressSettings,
) -> Result<PlaneStressState, MaterialError> {
    generic_state(material, in_plane, Some(initial_transverse), settings)
}

fn generic_state<L: FiniteStrainPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    initial_transverse: Option<[f64; 3]>,
    settings: PlaneStressSettings,
) -> Result<PlaneStressState, MaterialError> {
    plane_stress_state_with(
        in_plane,
        settings,
        initial_transverse,
        |f| material.stress_and_potential(f),
        |f| {
            let tangent = material.consistent_tangent(f)?;
            Ok(TRANSVERSE.map(|i| TRANSVERSE.map(|j| tangent[i][j])))
        },
    )
}

/// Consistent Schur-condensed response of any fixed finite-strain potential.
/// Reference-deformation pullbacks and inelastic return-map tangents are used
/// before transverse elimination, so their chain rules remain intact.
pub fn plane_stress<L: FiniteStrainPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    settings: PlaneStressSettings,
) -> Result<PlaneStressResponse, MaterialError> {
    let state = plane_stress_state(material, in_plane, settings)?;
    let tangent = material.consistent_tangent(&state.deformation_gradient)?;
    plane_stress_response_from(state, tangent)
}

/// Explicitly seeded finite-compliance equilibrium and its consistent tangent.
/// See [`plane_stress_state_with_seed`] for the seed's physical convention.
pub fn plane_stress_with_seed<L: FiniteStrainPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    initial_transverse: [f64; 3],
    settings: PlaneStressSettings,
) -> Result<PlaneStressResponse, MaterialError> {
    let state = plane_stress_state_with_seed(material, in_plane, initial_transverse, settings)?;
    let tangent = material.consistent_tangent(&state.deformation_gradient)?;
    plane_stress_response_from(state, tangent)
}

impl Hyperelastic {
    /// Relax all three transverse tractions using bounded, residual-decreasing
    /// Newton steps, without forming the full 9×9 constitutive Hessian.
    ///
    /// The initial transverse column is `[0,0,1]`. If its tangent is unstable,
    /// try a volume-preserving thickness; if that target leaves the material
    /// domain, bisect toward it from the admitted thickness to find a stable
    /// seed. At most `max_backtracks` candidate seeds are tried. Every trial
    /// retains positive orientation and the original principal-stretch domain;
    /// no stiffness or deformation is clipped. An unresolved stable seed or
    /// exhausted equilibrium budget refuses.
    /// Iteration and backtracking are bounded; spatial callers poll cancellation
    /// between quadrature points. No history or warm-start state is mutated.
    pub fn plane_stress_state(
        &self,
        in_plane: &[f64; 4],
        settings: PlaneStressSettings,
    ) -> Result<PlaneStressState, MaterialError> {
        plane_stress_state_with(
            in_plane,
            settings,
            None,
            |f| {
                self.check_state(f)?;
                let (energy, stress) = fs_ad::gradient::<9>(*f, |x| self.energy(&x));
                super::mixed::checked_stress(energy, stress)
            },
            |f| Ok(self.transverse_tangent(f)),
        )
    }

    /// Traction-free energy/stress and the exact Schur-complement tangent of
    /// the locally converged transverse equilibrium.
    pub fn plane_stress(
        &self,
        in_plane: &[f64; 4],
        settings: PlaneStressSettings,
    ) -> Result<PlaneStressResponse, MaterialError> {
        let state = self.plane_stress_state(in_plane, settings)?;
        let tangent = self.tangent(&state.deformation_gradient)?;
        plane_stress_response_from(state, tangent)
    }
}

// One bounded Newton/line-search owner. The native law supplies its compact
// three-direction AD block; opaque potential owners supply their exact tangent.
fn plane_stress_state_with(
    in_plane: &[f64; 4],
    settings: PlaneStressSettings,
    initial_transverse: Option<[f64; 3]>,
    point: impl Fn(&[f64; 9]) -> Result<([f64; 9], f64), MaterialError>,
    transverse: impl Fn(&[f64; 9]) -> Result<[[f64; 3]; 3], MaterialError>,
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
    if in_plane.iter().any(|v| !v.is_finite()) || !(determinant.is_finite() && determinant > 0.0) {
        return Err(state_error(
            "finite orientation-preserving in-plane gradient required",
        ));
    }
    let seed = initial_transverse.unwrap_or([0.0, 0.0, 1.0]);
    if seed.iter().any(|v| !v.is_finite()) || seed[2] <= 0.0 {
        return Err(state_error(
            "finite transverse seed with positive thickness required",
        ));
    }
    let mut f = [
        in_plane[0],
        in_plane[1],
        seed[0],
        in_plane[2],
        in_plane[3],
        seed[1],
        0.0,
        0.0,
        seed[2],
    ];
    let supplied = if initial_transverse.is_some() {
        match point(&f).and_then(|(stress, energy)| {
            transverse(&f)
                .and_then(PositiveBlock::new)
                .map(|block| (stress, energy, block))
        }) {
            Ok(data) => Some(data),
            Err(MaterialError::State { .. }) => {
                // A natural-volume guess is not an imposed constitutive mode.
                // For example, low-bulk NH can admit unit thickness while
                // det(G)/det(a) lies outside its elastic-stretch domain.
                f[2] = 0.0;
                f[5] = 0.0;
                f[8] = 1.0;
                None
            }
            Err(error) => return Err(error),
        }
    } else {
        None
    };
    let (mut p, mut energy, initial_block) = match supplied {
        Some((stress, energy, block)) => (stress, energy, Ok(block)),
        None => {
            let (stress, energy) = point(&f)?;
            (stress, energy, transverse(&f).and_then(PositiveBlock::new))
        }
    };
    let mut block = match initial_block {
        Ok(block) => block,
        Err(initial_error) => {
            // A poor seed need not be a physical instability. In particular,
            // log-squared bulk energy can have negative thickness curvature
            // at large J even though its traction-free root is stable.
            let target = 1.0 / determinant;
            if !(target.is_finite() && target > 0.0) || target == 1.0 {
                return Err(initial_error);
            }
            let (mut toward, mut admitted) = (target, 1.0);
            let mut candidate = target;
            let mut seed = None;
            for _ in 0..settings.max_backtracks {
                let mut trial = f;
                trial[8] = candidate;
                match point(&trial) {
                    Ok((stress, value)) => {
                        if let Ok(stable) = transverse(&trial).and_then(PositiveBlock::new) {
                            seed = Some((trial, stress, value, stable));
                            break;
                        }
                        // Admissible but unstable: continue toward J=1.
                        admitted = candidate;
                    }
                    Err(MaterialError::State { .. }) => toward = candidate,
                    Err(error) => return Err(error),
                }
                let next = 0.5 * toward + 0.5 * admitted;
                if next == candidate {
                    break;
                }
                candidate = next;
            }
            let Some((gradient, stress, value, stable)) = seed else {
                return Err(state_error(
                    "no admissible stable transverse seed within the initialization budget",
                ));
            };
            f = gradient;
            p = stress;
            energy = value;
            stable
        }
    };
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
                match point(&trial) {
                    Ok((stress, value)) => {
                        if relative_residual(&stress, scale) <= (1.0 - 1e-4 * alpha) * residual {
                            // Every accepted iterate remains on a stable
                            // transverse branch, not merely a traction root.
                            if let Ok(next) = transverse(&trial).and_then(PositiveBlock::new) {
                                accepted = Some((trial, stress, value, next));
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
        let Some((next_f, next_p, next_energy, next_block)) = accepted else {
            return Err(state_error(
                "transverse line search exhausted its budget or material domain",
            ));
        };
        f = next_f;
        p = next_p;
        energy = next_energy;
        block = next_block;
    }
    Err(state_error(
        "transverse Newton solve exhausted its iteration budget",
    ))
}

pub(crate) fn plane_stress_response_from(
    state: PlaneStressState,
    a: [[f64; 9]; 9],
) -> Result<PlaneStressResponse, MaterialError> {
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
