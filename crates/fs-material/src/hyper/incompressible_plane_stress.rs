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

use super::{
    Hyperelastic,
    cards::{HyperelasticMaterial, ResolvedHyperelastic},
    mixed::RemainderResponse,
    plane_stress::PlaneStressSettings,
    reference::ReferencedHyperelastic,
};
use crate::MaterialError;
use fs_ad::{Dual, Dual64, Real};

const IN_PLANE: [usize; 4] = [0, 1, 3, 4];
const SHEAR: [usize; 2] = [2, 5];

/// Converged exactly incompressible plane-stress energy and stress.
#[derive(Clone, Copy, Debug)]
pub struct IncompressiblePlaneStressState {
    /// Full row-major deformation gradient, with
    /// `F22=J_reference/det(F_in_plane)`.
    pub deformation_gradient: [f64; 9],
    /// Reduced remainder potential per original reference volume (Pa = J/m³).
    /// This is stored energy for elasticity. For a frozen inelastic owner it
    /// is an incremental potential, not physical storage or heat.
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
    /// Full row-major deformation gradient, with
    /// `F22=J_reference/det(F_in_plane)`.
    pub deformation_gradient: [f64; 9],
    /// Reduced remainder potential per original reference volume (Pa = J/m³).
    /// Frozen inelastic potentials are not physical stored energy or heat.
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

/// Exact-volume plane stress and its held-material kinematic derivatives.
///
/// Temperature, the natural reference (including its volume), constitutive
/// coefficients and accepted internal history are fixed. A frozen inelastic
/// owner supplies its consistent directional branch, not a complete material
/// time derivative. Independent source or history rates are not inferred.
#[derive(Clone, Copy, Debug)]
pub struct IncompressiblePlaneStressKinematicResponse {
    /// Original exact-volume response, including full pressure-corrected stress.
    pub response: IncompressiblePlaneStressResponse,
    /// d[F02,F12,F22]/d[F00,F01,F10,F11], including shear stationarity and
    /// the exact derivative of F22=J_reference/det(F_in_plane).
    pub transverse_deformation_tangent: [[f64; 4]; 3],
    /// Derivative of the native pressure multiplier with respect to the four
    /// in-plane deformation entries [Pa]. This retains the remainder's gauge.
    pub pressure_tangent_pa: [f64; 4],
}

impl IncompressiblePlaneStressKinematicResponse {
    /// Full Fdot from the supplied in-plane rate [1/s], preserving the held
    /// reference-volume constraint. The remaining bottom-row rates are zero.
    pub fn deformation_gradient_rate(
        &self,
        in_plane_rate: [f64; 4],
    ) -> Result<[f64; 9], MaterialError> {
        if in_plane_rate.iter().any(|v| !v.is_finite()) {
            return Err(state_error("nonfinite in-plane deformation rate"));
        }
        let mut rate = [0.0; 9];
        for col in 0..4 {
            rate[IN_PLANE[col]] = in_plane_rate[col];
        }
        for (row, index) in [2, 5, 8].into_iter().enumerate() {
            rate[index] = (0..4).fold(0.0, |sum, col| {
                self.transverse_deformation_tangent[row][col].mul_add(in_plane_rate[col], sum)
            });
        }
        if rate.iter().any(|v| !v.is_finite()) {
            return Err(state_error("unrepresentable exact-volume deformation rate"));
        }
        Ok(rate)
    }

    /// Rate of the compression-positive pressure multiplier [Pa/s] along the
    /// same held-material exact-volume branch. This is not a source-induced rate.
    pub fn pressure_rate(&self, in_plane_rate: [f64; 4]) -> Result<f64, MaterialError> {
        if in_plane_rate.iter().any(|v| !v.is_finite()) {
            return Err(state_error("nonfinite in-plane deformation rate"));
        }
        let rate = (0..4).fold(0.0, |sum, col| {
            self.pressure_tangent_pa[col].mul_add(in_plane_rate[col], sum)
        });
        if !rate.is_finite() {
            return Err(state_error("unrepresentable pressure rate"));
        }
        Ok(rate)
    }
}

fn state_error(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("incompressible hyperelastic plane stress: {what}"),
    }
}

// Six independent coordinates. The determinant constraint is differentiated
// through this map, not imposed after differentiating a compressible tangent.
fn constrained_gradient<T: Real>(x: &[T; 6], reference_volume: f64) -> [T; 9] {
    let thickness = T::from_f64(reference_volume) / (x[0] * x[3] - x[1] * x[2]);
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

/// A fixed-history, bulk-free potential with explicit elastic incompressibility.
///
/// This capability is separate from a finite-pressure Legendre transform.
/// Coupled Neo-Hookean has an admitted exact plane-stress remainder even though
/// the finite-compliance mixed formulation deliberately refuses that split.
/// Implementations retain their reference mapping and all accepted material
/// history. They must not evaluate a bulk penalty inside the remainder.
pub trait IncompressiblePlaneStressPotential {
    /// Required total volume ratio for incompressible elastic deformation.
    /// A native law returns one; a prescribed relaxed reference returns its
    /// determinant in the caller's original reference configuration.
    fn incompressible_reference_volume_ratio(&self) -> f64;

    /// Original-volume remainder stress and potential at fixed history.
    fn incompressible_remainder_stress_and_potential(
        &self,
        deformation: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError>;

    /// The same remainder and its consistent Hessian, including every
    /// reference-map chain rule and any frozen inelastic return map.
    fn incompressible_remainder_response(
        &self,
        deformation: &[f64; 9],
    ) -> Result<RemainderResponse, MaterialError>;
}

impl IncompressiblePlaneStressPotential for Hyperelastic {
    fn incompressible_reference_volume_ratio(&self) -> f64 {
        1.0
    }

    fn incompressible_remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        self.check_state(f)?;
        let (energy, stress) = fs_ad::gradient::<9>(*f, |x| self.isochoric_energy(&x));
        super::mixed::checked_stress(energy, stress)
    }

    fn incompressible_remainder_response(
        &self,
        f: &[f64; 9],
    ) -> Result<RemainderResponse, MaterialError> {
        self.check_state(f)?;
        super::mixed::response_from_energy(f, |x| self.isochoric_energy(x))
    }
}

impl IncompressiblePlaneStressPotential for ResolvedHyperelastic {
    fn incompressible_reference_volume_ratio(&self) -> f64 {
        self.law().incompressible_reference_volume_ratio()
    }

    fn incompressible_remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        self.law().incompressible_remainder_stress_and_potential(f)
    }

    fn incompressible_remainder_response(
        &self,
        f: &[f64; 9],
    ) -> Result<RemainderResponse, MaterialError> {
        self.law().incompressible_remainder_response(f)
    }
}

impl<L: HyperelasticMaterial + IncompressiblePlaneStressPotential + ?Sized>
    IncompressiblePlaneStressPotential for ReferencedHyperelastic<'_, L>
{
    fn incompressible_reference_volume_ratio(&self) -> f64 {
        self.reference_deformation().volume_ratio()
            * self.material().incompressible_reference_volume_ratio()
    }

    fn incompressible_remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        let reference = self.reference_deformation();
        let (stress, energy) = self
            .material()
            .incompressible_remainder_stress_and_potential(&reference.elastic_deformation(f))?;
        super::mixed::checked_stress(
            reference.volume_ratio() * energy,
            reference.pullback_piola(&stress),
        )
    }

    fn incompressible_remainder_response(
        &self,
        f: &[f64; 9],
    ) -> Result<RemainderResponse, MaterialError> {
        let reference = self.reference_deformation();
        let elastic = self
            .material()
            .incompressible_remainder_response(&reference.elastic_deformation(f))?;
        let (first_piola_pa, potential_pa) = super::mixed::checked_stress(
            reference.volume_ratio() * elastic.potential_pa,
            reference.pullback_piola(&elastic.first_piola_pa),
        )?;
        let tangent_pa = reference.pullback_tangent(&elastic.tangent_pa);
        if tangent_pa.iter().flatten().any(|x| !x.is_finite()) {
            return Err(state_error(
                "unrepresentable reference-pulled remainder tangent",
            ));
        }
        Ok(RemainderResponse {
            potential_pa,
            first_piola_pa,
            tangent_pa,
        })
    }
}

/// Relax transverse shear with `det(F)=J_reference`, preserving the supplied
/// potential's original reference map, constitutive domain and frozen history.
/// No bulk penalty is evaluated. Inelastic callers must separately commit
/// the accepted history and account for physical storage and dissipation.
pub fn incompressible_plane_stress_state<L: IncompressiblePlaneStressPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    settings: PlaneStressSettings,
) -> Result<IncompressiblePlaneStressState, MaterialError> {
    generic_state(material, in_plane, None, settings).map(|(state, _)| state)
}

/// Exact-volume relaxation from a retained or caller-supplied transverse shear
/// `[F02,F12]`. Thickness is always recomputed as `J_reference/det(a)`; a warm
/// start never changes the incompressibility constraint or material history.
/// Nonfinite seed entries refuse. A finite seed outside the material domain or
/// stable shear branch falls back to zero shear once, before Newton begins.
pub fn incompressible_plane_stress_state_with_seed<
    L: IncompressiblePlaneStressPotential + ?Sized,
>(
    material: &L,
    in_plane: &[f64; 4],
    initial_shear: [f64; 2],
    settings: PlaneStressSettings,
) -> Result<IncompressiblePlaneStressState, MaterialError> {
    generic_state(material, in_plane, Some(initial_shear), settings).map(|(state, _)| state)
}

/// Exact-volume plane stress with the consistent constrained-energy tangent.
/// The chain rule includes both derivatives of `F22=J_reference/det(a)` before
/// eliminating the two stationary transverse shears. This avoids subtracting
/// a large bulk penalty from a compressible tangent.
pub fn incompressible_plane_stress<L: IncompressiblePlaneStressPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    settings: PlaneStressSettings,
) -> Result<IncompressiblePlaneStressResponse, MaterialError> {
    let (state, remainder) = generic_state(material, in_plane, None, settings)?;
    generic_response(state, remainder, in_plane)
}

/// Seeded exact-volume equilibrium with the same consistent condensed tangent.
/// See [`incompressible_plane_stress_state_with_seed`] for seed admission.
pub fn incompressible_plane_stress_with_seed<L: IncompressiblePlaneStressPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    initial_shear: [f64; 2],
    settings: PlaneStressSettings,
) -> Result<IncompressiblePlaneStressResponse, MaterialError> {
    let (state, remainder) = generic_state(material, in_plane, Some(initial_shear), settings)?;
    generic_response(state, remainder, in_plane)
}

/// Exact-volume response with all relaxed transverse and pressure derivatives.
///
/// The optional seed contains only the two transverse shears; thickness always
/// obeys the exact volume constraint. No penalty law or evolving reference,
/// temperature, material coefficient or internal-history rate is evaluated.
pub fn incompressible_plane_stress_with_kinematics<
    L: IncompressiblePlaneStressPotential + ?Sized,
>(
    material: &L,
    in_plane: &[f64; 4],
    initial_shear: Option<[f64; 2]>,
    settings: PlaneStressSettings,
) -> Result<IncompressiblePlaneStressKinematicResponse, MaterialError> {
    let (state, remainder) = generic_state(material, in_plane, initial_shear, settings)?;
    let remainder = remainder.ok_or_else(|| state_error("missing accepted remainder tangent"))?;
    kinematic_response(state, remainder, in_plane)
}

pub(crate) fn generic_response(
    state: IncompressiblePlaneStressState,
    remainder: Option<RemainderResponse>,
    in_plane: &[f64; 4],
) -> Result<IncompressiblePlaneStressResponse, MaterialError> {
    let remainder = remainder.ok_or_else(|| state_error("missing accepted remainder tangent"))?;
    let (hessian, _) = constrained_hessian(&state, &remainder, in_plane);
    incompressible_response_from(state, hessian)
}

// Differentiate the exact constraint before condensing the relaxed shears.
// Retaining the thickness gradient also supplies the kinematic closure.
fn constrained_hessian(
    state: &IncompressiblePlaneStressState,
    remainder: &RemainderResponse,
    in_plane: &[f64; 4],
) -> ([[f64; 6]; 6], [f64; 4]) {
    let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
    let thickness = state.deformation_gradient[8];
    let determinant_gradient = [in_plane[3], -in_plane[2], -in_plane[1], in_plane[0]];
    let thickness_gradient = determinant_gradient.map(|v| -(thickness / determinant) * v);
    let determinant_hessian = [
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0, 0.0],
        [0.0, -1.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
    ];
    let independent = [0, 1, 3, 4, 2, 5];
    let mut hessian = [[0.0; 6]; 6];
    for row in 0..6 {
        let dt_row = if row < 4 {
            thickness_gradient[row]
        } else {
            0.0
        };
        for col in 0..6 {
            let dt_col = if col < 4 {
                thickness_gradient[col]
            } else {
                0.0
            };
            let d2t = if row < 4 && col < 4 {
                -(dt_row * determinant_gradient[col]
                    + dt_col * determinant_gradient[row]
                    + thickness * determinant_hessian[row][col])
                    / determinant
            } else {
                0.0
            };
            let (i, j) = (independent[row], independent[col]);
            hessian[row][col] = remainder.tangent_pa[i][j]
                + remainder.tangent_pa[i][8] * dt_col
                + remainder.tangent_pa[8][j] * dt_row
                + remainder.tangent_pa[8][8] * dt_row * dt_col
                + remainder.first_piola_pa[8] * d2t;
        }
    }
    (hessian, thickness_gradient)
}

pub(crate) fn kinematic_response(
    state: IncompressiblePlaneStressState,
    remainder: RemainderResponse,
    in_plane: &[f64; 4],
) -> Result<IncompressiblePlaneStressKinematicResponse, MaterialError> {
    let (hessian, thickness_gradient) = constrained_hessian(&state, &remainder, in_plane);
    let (response, shear_gradient) = condensed_response(state, hessian)?;
    let transverse_deformation_tangent = [shear_gradient[0], shear_gradient[1], thickness_gradient];
    let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
    let determinant_gradient = [in_plane[3], -in_plane[2], -in_plane[1], in_plane[0]];
    let pressure_tangent_pa = core::array::from_fn(|col| {
        // p=Pr22/d, so pdot=(A22:Fdot-p*ddot)/d. The held reference
        // volume has zero rate; no artificial bulk stiffness enters.
        let remainder_rate = [2, 5, 8].into_iter().enumerate().fold(
            remainder.tangent_pa[8][IN_PLANE[col]],
            |sum, (row, index)| {
                remainder.tangent_pa[8][index]
                    .mul_add(transverse_deformation_tangent[row][col], sum)
            },
        );
        (remainder_rate - state.pressure_pa * determinant_gradient[col]) / determinant
    });
    if transverse_deformation_tangent
        .iter()
        .flatten()
        .any(|v| !v.is_finite())
        || pressure_tangent_pa.iter().any(|v| !v.is_finite())
    {
        return Err(state_error("unrepresentable exact-volume motion tangent"));
    }
    Ok(IncompressiblePlaneStressKinematicResponse {
        response,
        transverse_deformation_tangent,
        pressure_tangent_pa,
    })
}

pub(crate) fn generic_state<L: IncompressiblePlaneStressPotential + ?Sized>(
    material: &L,
    in_plane: &[f64; 4],
    initial_shear: Option<[f64; 2]>,
    settings: PlaneStressSettings,
) -> Result<(IncompressiblePlaneStressState, Option<RemainderResponse>), MaterialError> {
    incompressible_state_with(
        in_plane,
        settings,
        material.incompressible_reference_volume_ratio(),
        initial_shear,
        |f| {
            let response = material.incompressible_remainder_response(f)?;
            if !response.potential_pa.is_finite()
                || response.first_piola_pa.iter().any(|v| !v.is_finite())
            {
                return Err(state_error("nonfinite remainder potential or stress"));
            }
            Ok(ShearEvaluation {
                stress: SHEAR.map(|i| response.first_piola_pa[i]),
                block: PositiveShear::new(SHEAR.map(|i| SHEAR.map(|j| response.tangent_pa[i][j])))?,
                response: Some(response),
            })
        },
        |f| material.incompressible_remainder_stress_and_potential(f),
    )
}

// Scalar-normalized Cholesky of the two eliminated shear coordinates.
// An unstable block refuses; no regularization or stiffness floor is used.
pub(crate) struct PositiveShear {
    lower: [[f64; 2]; 2],
    scale: f64,
}

// Generic inelastic owners retain the prepared point through acceptance so
// condensation does not repeat their local history solve. The native elastic
// fast path forms only a two-direction shear block during local Newton.
struct ShearEvaluation {
    stress: [f64; 2],
    block: PositiveShear,
    response: Option<RemainderResponse>,
}

impl PositiveShear {
    pub(crate) fn new(a: [[f64; 2]; 2]) -> Result<Self, MaterialError> {
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

    pub(crate) fn solve(&self, rhs: [f64; 2]) -> Result<[f64; 2], MaterialError> {
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
        incompressible_state_with(
            in_plane,
            settings,
            1.0,
            None,
            |f| {
                let (stress, block) = self.incompressible_shear_data(f)?;
                Ok(ShearEvaluation {
                    stress,
                    block,
                    response: None,
                })
            },
            |f| self.incompressible_remainder_stress_and_potential(f),
        )
        .map(|(state, _)| state)
    }
}

fn incompressible_state_with(
    in_plane: &[f64; 4],
    settings: PlaneStressSettings,
    reference_volume: f64,
    initial_shear: Option<[f64; 2]>,
    shear: impl Fn(&[f64; 9]) -> Result<ShearEvaluation, MaterialError>,
    point: impl Fn(&[f64; 9]) -> Result<([f64; 9], f64), MaterialError>,
) -> Result<(IncompressiblePlaneStressState, Option<RemainderResponse>), MaterialError> {
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
    if !(reference_volume.is_finite() && reference_volume > 0.0) {
        return Err(state_error(
            "finite positive incompressible reference volume required",
        ));
    }
    let determinant = in_plane[0] * in_plane[3] - in_plane[1] * in_plane[2];
    if in_plane.iter().any(|v| !v.is_finite()) || !(determinant.is_finite() && determinant > 0.0) {
        return Err(state_error(
            "finite orientation-preserving in-plane gradient required",
        ));
    }
    let seed = initial_shear.unwrap_or([0.0, 0.0]);
    if seed.iter().any(|value| !value.is_finite()) {
        return Err(state_error("finite transverse shear seed required"));
    }
    let mut f = constrained_gradient(
        &[
            in_plane[0],
            in_plane[1],
            in_plane[2],
            in_plane[3],
            seed[0],
            seed[1],
        ],
        reference_volume,
    );
    let mut data = match shear(&f) {
        Ok(data) => data,
        Err(MaterialError::State { .. }) if initial_shear.is_some() && seed != [0.0, 0.0] => {
            f[2] = 0.0;
            f[5] = 0.0;
            shear(&f)?
        }
        Err(error) => return Err(error),
    };
    let gradient_norm = in_plane[0]
        .hypot(in_plane[1])
        .hypot(in_plane[2])
        .hypot(in_plane[3]);
    let scale = data.block.traction_scale()? * gradient_norm.max(1.0);
    if !(scale.is_finite() && scale > 0.0) {
        return Err(state_error("unrepresentable traction normalization"));
    }
    for iteration in 0..=settings.max_iterations {
        let residual = (data.stress[0] / scale).hypot(data.stress[1] / scale);
        if residual <= settings.relative_tolerance {
            let (remainder_piola, energy) = match data.response {
                Some(response) => (response.first_piola_pa, response.potential_pa),
                None => point(&f)?,
            };
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
            let full_piola = core::array::from_fn(|i| remainder_piola[i] - pressure * cofactor[i]);
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
            return Ok((
                IncompressiblePlaneStressState {
                    deformation_gradient: f,
                    energy_density_pa: energy,
                    piola: IN_PLANE.map(|i| full_piola[i]),
                    full_piola,
                    pressure_pa: pressure,
                    relative_transverse_residual: final_residual,
                    iterations: iteration,
                },
                data.response,
            ));
        }
        if iteration == settings.max_iterations || !residual.is_finite() {
            break;
        }
        let direction = data.block.solve(data.stress.map(|v| -v))?;
        let mut alpha = 1.0;
        let mut accepted = None;
        for _ in 0..settings.max_backtracks {
            let mut trial = f;
            for i in 0..2 {
                trial[SHEAR[i]] += alpha * direction[i];
            }
            match shear(&trial) {
                Ok(next) => {
                    let next_residual = (next.stress[0] / scale).hypot(next.stress[1] / scale);
                    if next_residual <= (1.0 - 1e-4 * alpha) * residual {
                        accepted = Some((trial, next));
                        break;
                    }
                }
                Err(MaterialError::State { .. }) => (),
                Err(error) => return Err(error),
            }
            alpha *= 0.5;
        }
        let Some((next_f, next)) = accepted else {
            return Err(state_error(
                "transverse shear line search exhausted its budget or material domain",
            ));
        };
        f = next_f;
        data = next;
    }
    Err(state_error(
        "transverse shear Newton solve exhausted its iteration budget",
    ))
}

impl Hyperelastic {
    /// Exact-volume kinematics and pressure derivative on the stable relaxed
    /// branch, at fixed temperature, reference and accepted material history.
    /// Native Newton uses its two-direction shear block; the full remainder
    /// tangent is formed once at acceptance for the constrained derivatives.
    pub fn incompressible_plane_stress_with_kinematics(
        &self,
        in_plane: &[f64; 4],
        settings: PlaneStressSettings,
    ) -> Result<IncompressiblePlaneStressKinematicResponse, MaterialError> {
        let state = self.incompressible_plane_stress_state(in_plane, settings)?;
        let remainder = self.incompressible_remainder_response(&state.deformation_gradient)?;
        kinematic_response(state, remainder, in_plane)
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
            let w = self.isochoric_energy(&constrained_gradient(&vars, 1.0));
            for row in 0..6 {
                hessian[row][col] = w.eps[0].eps[row];
            }
        }
        incompressible_response_from(state, hessian)
    }
}

fn incompressible_response_from(
    state: IncompressiblePlaneStressState,
    hessian: [[f64; 6]; 6],
) -> Result<IncompressiblePlaneStressResponse, MaterialError> {
    condensed_response(state, hessian).map(|(response, _)| response)
}

fn condensed_response(
    state: IncompressiblePlaneStressState,
    hessian: [[f64; 6]; 6],
) -> Result<(IncompressiblePlaneStressResponse, [[f64; 4]; 2]), MaterialError> {
    if hessian.iter().flatten().any(|v| !v.is_finite()) {
        return Err(state_error("nonfinite constrained-energy Hessian"));
    }
    let block = PositiveShear::new([
        [hessian[4][4], hessian[4][5]],
        [hessian[5][4], hessian[5][5]],
    ])?;
    let mut tangent = [[0.0; 4]; 4];
    let mut shear_gradient = [[0.0; 4]; 2];
    for col in 0..4 {
        let extension = block.solve([hessian[4][col], hessian[5][col]])?;
        for row in 0..2 {
            shear_gradient[row][col] = -extension[row];
        }
        for row in 0..4 {
            tangent[row][col] =
                hessian[row][col] - hessian[row][4] * extension[0] - hessian[row][5] * extension[1];
        }
    }
    if tangent.iter().flatten().any(|v| !v.is_finite()) {
        return Err(state_error("unrepresentable condensed tangent"));
    }
    Ok((
        IncompressiblePlaneStressResponse {
            deformation_gradient: state.deformation_gradient,
            energy_density_pa: state.energy_density_pa,
            piola: state.piola,
            full_piola: state.full_piola,
            pressure_pa: state.pressure_pa,
            tangent,
            relative_transverse_residual: state.relative_transverse_residual,
            iterations: state.iterations,
        },
        shear_gradient,
    ))
}
