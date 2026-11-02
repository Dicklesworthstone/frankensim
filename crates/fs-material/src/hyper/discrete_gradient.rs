//! Objective two-point elastic work and its consistent endpoint derivative.
//!
//! The Gonzalez discrete derivative is formed in the right Cauchy–Green
//! strain, not in an interpolated deformation. With `C_m=(C_0+C_1)/2`,
//! `d=C_1-C_0`, and `g=dW/dC`, it is
//! `g_d=g(C_m)+(W(C_1)-W(C_0)-g(C_m):d)*d/(d:d)`.
//! The algorithmic Piola stress is `P_d=(F_0+F_1)*g_d`. Its work equals
//! the endpoint energy change, while `P_d*F_m^T` is symmetric. These are
//! the local energy and angular-momentum identities used by a conserving
//! spatial/time discretization; this material owner is not a time integrator.
//!
//! Mandel coordinates preserve the Frobenius metric. A Cholesky representative
//! of C feeds the ORIGINAL generic-scalar energy and native AD, including
//! anisotropic fibres and repeated-stretch Ogden response. The consistent
//! endpoint Jacobian is generally nonsymmetric and must not be treated as
//! an elastic Hessian. Only first/second native derivatives are required.
//!
//! At an explicitly declared small strain increment, the uncorrected strain
//! midpoint avoids dividing cancellation by a vanishing norm. Its actual work
//! defect is reported and admitted against explicit absolute/relative bounds;
//! a silently clamped exact-energy claim is never made. The same mechanism
//! supplies an objective discrete derivative of J for mixed constraints.
//!
//! Reference: O. Gonzalez, CMAME 190 (2000), 1763–1783,
//! doi:10.1016/S0045-7825(00)00189-4, equations (10)–(12).

use fs_ad::{Dual, Dual64, Real};

use super::mixed::{MixedFiniteStrainPotential, cofactor};
use super::{Hyperelastic, det3};
use crate::MaterialError;

/// Original finite-bulk storage and its work-conjugate interval pressure.
pub mod bulk;

const SQRT_2: f64 = core::f64::consts::SQRT_2;

/// The complete original elastic energy or its admitted bulk-free remainder.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EnergyPart {
    Full,
    /// Native mixed split; coupled Neo-Hookean is not admitted here.
    Remainder,
}

/// Explicit cancellation resolution and local work admission. No defaults.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscreteGradientSettings {
    /// Use the strain-midpoint limit when `|C1-C0| <= this*max(|C0|,|C1|)`.
    pub minimum_relative_strain_increment: f64,
    /// Allowed absolute constitutive work defect per reference volume [Pa].
    pub absolute_work_tolerance_pa: f64,
    /// Additional allowance times max(|W0|, |W1|, |algorithmic work|).
    pub relative_work_tolerance: f64,
}

impl DiscreteGradientSettings {
    fn validate(self) -> Result<(), MaterialError> {
        if ![
            self.minimum_relative_strain_increment,
            self.absolute_work_tolerance_pa,
            self.relative_work_tolerance,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
            || self.minimum_relative_strain_increment > 1.0
        {
            return Err(MaterialError::Parameters {
                what: "finite nonnegative discrete-gradient controls and strain threshold <= 1 required"
                    .into(),
            });
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscreteGradientState {
    pub energy_before_pa: f64,
    pub energy_after_pa: f64,
    /// Algorithmic force for this pair of states, not endpoint elastic stress.
    pub first_piola_pa: [f64; 9],
    pub work_pa: f64,
    pub work_defect_pa: f64,
    pub used_midpoint_limit: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscreteGradientResponse {
    pub state: DiscreteGradientState,
    /// dP_algorithmic/dF1 with F0 fixed. In general NONSYMMETRIC.
    pub endpoint_tangent_pa: [[f64; 9]; 9],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscreteVolumeState {
    pub volume_ratio_before: f64,
    pub volume_ratio_after: f64,
    /// Work-conjugate volume derivative: Gd:(F1-F0) = J1-J0 to its budget.
    pub discrete_cofactor: [f64; 9],
    /// Derivative of endpoint J1. Distinct from Gd at finite increments.
    pub endpoint_cofactor: [f64; 9],
    pub volume_work: f64,
    pub volume_work_defect: f64,
    pub used_midpoint_limit: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiscreteVolumeResponse {
    pub state: DiscreteVolumeState,
    /// Derivative of discrete_cofactor with respect to F1, with F0 fixed.
    pub endpoint_jacobian: [[f64; 9]; 9],
}

/// Immutable accepted deformation; repeated Newton trials cannot advance it.
#[derive(Clone, Copy, Debug)]
pub struct DiscreteGradientPoint<'a> {
    law: &'a Hyperelastic,
    part: EnergyPart,
    accepted_f: [f64; 9],
    accepted_c: [f64; 6],
    energy_before_pa: f64,
    volume_ratio_before: f64,
    settings: DiscreteGradientSettings,
}

#[derive(Clone, Copy)]
enum Potential<'a> {
    Material(&'a Hyperelastic, EnergyPart),
    Volume,
}

fn error(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("objective discrete gradient: {what}"),
    }
}

fn mandel_c(f: &[f64; 9]) -> [f64; 6] {
    let c = |a: usize, b: usize| (0..3).fold(0.0, |s, i| f[3 * i + a].mul_add(f[3 * i + b], s));
    [
        c(0, 0),
        c(1, 1),
        c(2, 2),
        SQRT_2 * c(0, 1),
        SQRT_2 * c(0, 2),
        SQRT_2 * c(1, 2),
    ]
}

// Upper Cholesky R satisfies R^T R=C. No inverse of F_mid is used: that
// midpoint may be singular even when both endpoint deformations are valid.
fn representative<T: Real>(c: &[T; 6]) -> [T; 9] {
    let off = T::from_f64(1.0 / SQRT_2);
    let r00 = c[0].sqrt();
    let r01 = c[3] * off / r00;
    let r02 = c[4] * off / r00;
    let r11 = (c[1] - r01 * r01).sqrt();
    let r12 = (c[5] * off - r01 * r02) / r11;
    let r22 = (c[2] - r02 * r02 - r12 * r12).sqrt();
    [
        r00,
        r01,
        r02,
        T::zero(),
        r11,
        r12,
        T::zero(),
        T::zero(),
        r22,
    ]
}

impl Potential<'_> {
    fn value<T: Real>(self, c: &[T; 6]) -> T {
        let f = representative(c);
        match self {
            Self::Material(law, EnergyPart::Full) => law.energy(&f),
            Self::Material(law, EnergyPart::Remainder) => law.isochoric_energy(&f),
            Self::Volume => det3(&f),
        }
    }

    fn derivatives(
        self,
        c: [f64; 6],
        tangent: bool,
    ) -> Result<([f64; 6], [[f64; 6]; 6]), MaterialError> {
        let mut h = [[0.0; 6]; 6];
        let g = if tangent {
            let mut g = [0.0; 6];
            for col in 0..6 {
                let vars = core::array::from_fn(|row| Dual {
                    re: Dual64::variable(c[row], row),
                    eps: [Dual64::constant(if row == col { 1.0 } else { 0.0 })],
                });
                let w = self.value(&vars);
                if col == 0 {
                    g = w.re.eps;
                }
                for row in 0..6 {
                    h[row][col] = w.eps[0].eps[row];
                }
            }
            g
        } else {
            fs_ad::gradient::<6>(c, |x| self.value(&x)).1
        };
        if g.iter().chain(h.iter().flatten()).any(|v| !v.is_finite()) {
            return Err(error("nonfinite native strain derivative"));
        }
        Ok((g, h))
    }
}

fn second_piola(g: &[f64; 6]) -> [[f64; 3]; 3] {
    [
        [2.0 * g[0], SQRT_2 * g[3], SQRT_2 * g[4]],
        [SQRT_2 * g[3], 2.0 * g[1], SQRT_2 * g[5]],
        [SQRT_2 * g[4], SQRT_2 * g[5], 2.0 * g[2]],
    ]
}

fn c_direction(f: &[f64; 9], column: usize) -> [f64; 6] {
    let (i, a) = (column / 3, column % 3);
    let d = |b, c| {
        (if a == b { f[3 * i + c] } else { 0.0 }) + (if a == c { f[3 * i + b] } else { 0.0 })
    };
    [
        d(0, 0),
        d(1, 1),
        d(2, 2),
        SQRT_2 * d(0, 1),
        SQRT_2 * d(0, 2),
        SQRT_2 * d(1, 2),
    ]
}

struct PairResponse {
    value_after: f64,
    gradient: [f64; 9],
    tangent: [[f64; 9]; 9],
    work: f64,
    defect: f64,
    limit: bool,
}

impl<'a> DiscreteGradientPoint<'a> {
    pub fn try_new(
        law: &'a Hyperelastic,
        part: EnergyPart,
        accepted_f: [f64; 9],
        settings: DiscreteGradientSettings,
    ) -> Result<Self, MaterialError> {
        settings.validate()?;
        match part {
            EnergyPart::Full => {
                Hyperelastic::new(law.model.clone(), law.stretch_limit)?;
            }
            EnergyPart::Remainder => law.validate_mixed()?,
        }
        law.check_state(&accepted_f)?;
        let accepted_c = mandel_c(&accepted_f);
        let energy_before_pa = Potential::Material(law, part).value(&accepted_c);
        let volume_ratio_before = Potential::Volume.value(&accepted_c);
        if !energy_before_pa.is_finite()
            || !volume_ratio_before.is_finite()
            || volume_ratio_before <= 0.0
        {
            return Err(error("accepted strain energy or volume is unrepresentable"));
        }
        Ok(Self {
            law,
            part,
            accepted_f,
            accepted_c,
            energy_before_pa,
            volume_ratio_before,
            settings,
        })
    }

    pub fn accepted_deformation(&self) -> &[f64; 9] {
        &self.accepted_f
    }
    pub fn energy_before_pa(&self) -> f64 {
        self.energy_before_pa
    }

    pub fn state(&self, next_f: &[f64; 9]) -> Result<DiscreteGradientState, MaterialError> {
        let pair = self.evaluate(
            next_f,
            Potential::Material(self.law, self.part),
            self.energy_before_pa,
            self.settings.absolute_work_tolerance_pa,
            false,
        )?;
        Ok(self.material_state(pair))
    }

    pub fn response(&self, next_f: &[f64; 9]) -> Result<DiscreteGradientResponse, MaterialError> {
        let pair = self.evaluate(
            next_f,
            Potential::Material(self.law, self.part),
            self.energy_before_pa,
            self.settings.absolute_work_tolerance_pa,
            true,
        )?;
        let endpoint_tangent_pa = pair.tangent;
        Ok(DiscreteGradientResponse {
            state: self.material_state(pair),
            endpoint_tangent_pa,
        })
    }

    fn material_state(&self, pair: PairResponse) -> DiscreteGradientState {
        DiscreteGradientState {
            energy_before_pa: self.energy_before_pa,
            energy_after_pa: pair.value_after,
            first_piola_pa: pair.gradient,
            work_pa: pair.work,
            work_defect_pa: pair.defect,
            used_midpoint_limit: pair.limit,
        }
    }

    /// Dimensionless absolute work allowance for J; material Pa allowance is not used.
    pub fn volume_state(
        &self,
        next_f: &[f64; 9],
        absolute_volume_work_tolerance: f64,
    ) -> Result<DiscreteVolumeState, MaterialError> {
        let pair = self.evaluate(
            next_f,
            Potential::Volume,
            self.volume_ratio_before,
            absolute_volume_work_tolerance,
            false,
        )?;
        Ok(self.volume_result(next_f, pair))
    }

    pub fn volume_response(
        &self,
        next_f: &[f64; 9],
        absolute_volume_work_tolerance: f64,
    ) -> Result<DiscreteVolumeResponse, MaterialError> {
        let pair = self.evaluate(
            next_f,
            Potential::Volume,
            self.volume_ratio_before,
            absolute_volume_work_tolerance,
            true,
        )?;
        let endpoint_jacobian = pair.tangent;
        Ok(DiscreteVolumeResponse {
            state: self.volume_result(next_f, pair),
            endpoint_jacobian,
        })
    }

    fn volume_result(&self, next_f: &[f64; 9], pair: PairResponse) -> DiscreteVolumeState {
        DiscreteVolumeState {
            volume_ratio_before: self.volume_ratio_before,
            volume_ratio_after: pair.value_after,
            discrete_cofactor: pair.gradient,
            endpoint_cofactor: cofactor(next_f),
            volume_work: pair.work,
            volume_work_defect: pair.defect,
            used_midpoint_limit: pair.limit,
        }
    }

    fn evaluate(
        &self,
        next_f: &[f64; 9],
        potential: Potential<'_>,
        before: f64,
        absolute_tolerance: f64,
        with_tangent: bool,
    ) -> Result<PairResponse, MaterialError> {
        if !absolute_tolerance.is_finite() || absolute_tolerance < 0.0 {
            return Err(error("finite nonnegative absolute work allowance required"));
        }
        self.law.check_state(next_f)?;
        let next_c = mandel_c(next_f);
        let midpoint_c = core::array::from_fn(|i| 0.5 * self.accepted_c[i] + 0.5 * next_c[i]);
        let delta: [f64; 6] = core::array::from_fn(|i| next_c[i] - self.accepted_c[i]);
        let norm = |x: &[f64; 6]| x.iter().fold(0.0_f64, |s, x| s.hypot(*x));
        let dn = norm(&delta);
        let scale = norm(&self.accepted_c).max(norm(&next_c));
        let after = potential.value(&next_c);
        let increment = after - before;
        if !dn.is_finite() || !scale.is_finite() || !after.is_finite() || !increment.is_finite() {
            return Err(error(
                "endpoint strain, potential or increment is unrepresentable",
            ));
        }
        let limit = dn == 0.0 || dn / scale <= self.settings.minimum_relative_strain_increment;
        let (mid_g, mid_h) = potential.derivatives(midpoint_c, with_tangent)?;
        let mut g = mid_g;
        let mut dg = mid_h.map(|row| row.map(|v| 0.5 * v));
        if !limit {
            // Unit direction avoids forming |d|², which could under/overflow.
            let n = delta.map(|v| v / dn);
            let directional = (0..6).fold(0.0, |s, i| mid_g[i].mul_add(n[i], s));
            let correction = increment / dn - directional;
            for i in 0..6 {
                g[i] = correction.mul_add(n[i], mid_g[i]);
            }
            if with_tangent {
                let next_g = potential.derivatives(next_c, false)?.0;
                let alpha = correction / dn;
                for col in 0..6 {
                    let h_n = (0..6).fold(0.0, |s, k| mid_h[k][col].mul_add(n[k], s));
                    let derivative =
                        (next_g[col] - mid_g[col]) / dn - 0.5 * h_n - 2.0 * alpha * n[col];
                    for row in 0..6 {
                        dg[row][col] += n[row] * derivative + if row == col { alpha } else { 0.0 };
                    }
                }
            }
        }
        let fm: [f64; 9] = core::array::from_fn(|i| 0.5 * self.accepted_f[i] + 0.5 * next_f[i]);
        let s = second_piola(&g);
        let gradient: [f64; 9] = core::array::from_fn(|ij| {
            (0..3).fold(0.0, |sum, a| {
                fm[3 * (ij / 3) + a].mul_add(s[a][ij % 3], sum)
            })
        });
        let mut tangent = [[0.0; 9]; 9];
        if with_tangent {
            for col in 0..9 {
                let dc = c_direction(next_f, col);
                let ds = second_piola(&core::array::from_fn(|row| {
                    (0..6).fold(0.0, |sum, j| dg[row][j].mul_add(dc[j], sum))
                }));
                for row in 0..9 {
                    let (i, j) = (row / 3, row % 3);
                    tangent[row][col] = (if i == col / 3 {
                        0.5 * s[col % 3][j]
                    } else {
                        0.0
                    }) + (0..3)
                        .fold(0.0, |sum, a| fm[3 * i + a].mul_add(ds[a][j], sum));
                }
            }
        }
        let work = (0..9).fold(0.0, |sum, i| {
            gradient[i].mul_add(next_f[i] - self.accepted_f[i], sum)
        });
        let defect = work - increment;
        let allowance = absolute_tolerance
            + self.settings.relative_work_tolerance * before.abs().max(after.abs()).max(work.abs());
        if gradient
            .iter()
            .chain(tangent.iter().flatten())
            .any(|v| !v.is_finite())
            || !work.is_finite()
            || !defect.is_finite()
            || !allowance.is_finite()
        {
            return Err(error(
                "algorithmic stress, tangent or work is unrepresentable",
            ));
        }
        if defect.abs() > allowance {
            return Err(MaterialError::State {
                what: format!(
                    "objective discrete gradient work defect {defect:e} exceeds {allowance:e}"
                ),
            });
        }
        Ok(PairResponse {
            value_after: after,
            gradient,
            tangent,
            work,
            defect,
            limit,
        })
    }
}
