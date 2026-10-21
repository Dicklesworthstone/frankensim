//! Active fibre stress (plan §4 M0.3-W3): muscle activation as a PRESCRIBED stress along a reference fibre
//! direction, added to a passive hyperelastic law (M0.3-W1a).
//!
//! The active Cauchy stress is `sigma_a = T_max a f_L(lambda_f) m (x) m`, with `lambda_f = |F a0|`, the current
//! fibre direction `m = F a0 / lambda_f`, the activation `a` in `[0, 1]` and the isometric force-length factor
//! `f_L`. Its first Piola image is `P_a = J sigma_a F^-T = (J / lambda_f^2) T_max a f_L (F a0) (x) a0`.
//!
//! This is a PRESCRIBED, non-conservative input, not a stored energy: the passive law keeps its energy and
//! tangent symmetry, the active part adds an asymmetric tangent, and its work is booked as an external (player)
//! power input, `P_a : Fdot`, never as stored energy. Activation is a performance control; nothing in this law
//! sets it. Only the isometric law is implemented: there is no force-velocity factor (stated in the CONTRACT).
//!
//! Zero activation returns the passive law bitwise. A fibre stretch outside the declared force-length domain
//! refuses (`MaterialError::State`), as does an activation outside `[0, 1]`; at zero activation the domain is
//! not consulted because no active stress is evaluated.

use crate::hyper::{Hyperelastic, HyperelasticModel};
use crate::{MaterialAdmissibility, MaterialError};
use fs_ad::{Dual64, Real};

/// Largest number of knots of a tabulated force-length curve.
pub const MAX_FORCE_LENGTH_KNOTS: usize = 64;

/// Dimensionless isometric force-length factor `f_L(lambda_f)` with values in `[0, 1]`.
#[derive(Debug, Clone, PartialEq)]
pub enum ForceLength {
    /// `exp(-((l - stretch_opt) / width)^2)`.
    Gaussian {
        /// Fibre stretch of maximal force (positive).
        stretch_opt: f64,
        /// Width of the bell (positive).
        width: f64,
    },
    /// `max(0, 1 - ((l - stretch_opt) / half_width)^2)`.
    Parabola {
        /// Fibre stretch of maximal force (positive).
        stretch_opt: f64,
        /// Half width at which the force reaches zero (positive).
        half_width: f64,
    },
    /// Piecewise-linear through `(stretch, factor)` knots with strictly ascending stretches and factors in
    /// `[0, 1]`; the end factors extend constantly (the declared domain bounds the use).
    Table {
        /// The knots, at least two and at most [`MAX_FORCE_LENGTH_KNOTS`].
        knots: Vec<(f64, f64)>,
    },
}

impl ForceLength {
    /// The factor at a (generic-scalar) fibre stretch. Branches use the primal value, so a kink gives a
    /// one-sided derivative.
    pub fn factor<T: Real>(&self, l: T) -> T {
        match self {
            Self::Gaussian { stretch_opt, width } => {
                let x = (l - T::from_f64(*stretch_opt)) * T::from_f64(1.0 / *width);
                (T::zero() - x * x).exp()
            }
            Self::Parabola {
                stretch_opt,
                half_width,
            } => {
                let x = (l - T::from_f64(*stretch_opt)) * T::from_f64(1.0 / *half_width);
                let value = T::one() - x * x;
                if value.value() > 0.0 { value } else { T::zero() }
            }
            Self::Table { knots } => {
                if knots.is_empty() {
                    return T::zero();
                }
                let lv = l.value();
                let last = knots.len() - 1;
                if lv <= knots[0].0 {
                    return T::from_f64(knots[0].1);
                }
                if lv >= knots[last].0 {
                    return T::from_f64(knots[last].1);
                }
                let mut i = 0;
                while knots[i + 1].0 <= lv {
                    i += 1;
                }
                let (x0, f0) = knots[i];
                let (x1, f1) = knots[i + 1];
                T::from_f64(f0)
                    + T::from_f64(f1 - f0) * (l - T::from_f64(x0)) * T::from_f64(1.0 / (x1 - x0))
            }
        }
    }

    fn validate(&self) -> Result<(), MaterialError> {
        let bad = |what: &str| Err(MaterialError::Parameters { what: what.to_string() });
        match self {
            Self::Gaussian { stretch_opt, width } => {
                if !(stretch_opt.is_finite() && *stretch_opt > 0.0 && width.is_finite() && *width > 0.0)
                {
                    return bad("Gaussian force-length needs finite positive stretch_opt and width");
                }
            }
            Self::Parabola {
                stretch_opt,
                half_width,
            } => {
                if !(stretch_opt.is_finite()
                    && *stretch_opt > 0.0
                    && half_width.is_finite()
                    && *half_width > 0.0)
                {
                    return bad("parabolic force-length needs finite positive stretch_opt and half_width");
                }
            }
            Self::Table { knots } => {
                if knots.len() < 2 || knots.len() > MAX_FORCE_LENGTH_KNOTS {
                    return bad("tabulated force-length needs 2..=MAX_FORCE_LENGTH_KNOTS knots");
                }
                if !knots
                    .iter()
                    .all(|(l, f)| l.is_finite() && *l > 0.0 && f.is_finite() && (0.0..=1.0).contains(f))
                {
                    return bad("force-length knots need finite positive stretches and factors in [0, 1]");
                }
                if !knots.windows(2).all(|w| w[0].0 < w[1].0) {
                    return bad("force-length knot stretches must be strictly ascending");
                }
            }
        }
        Ok(())
    }
}

/// The active fibre stress law: the maximal isometric stress, the force-length family and the reference fibre.
#[derive(Debug, Clone, PartialEq)]
pub struct ActiveFibreStress {
    /// Maximal isometric active Cauchy stress [Pa] at full activation and `f_L = 1`.
    pub t_max_pa: f64,
    /// The isometric force-length factor.
    pub force_length: ForceLength,
    /// The reference fibre direction (unit vector, reference frame).
    pub fibre: [f64; 3],
    /// Inclusive fibre-stretch domain `(lo, hi)` of the force-length law; outside it the law refuses.
    pub stretch_domain: (f64, f64),
}

impl ActiveFibreStress {
    /// Construct with parameter checks.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] for a non-positive or non-finite `t_max_pa`, a non-unit or non-finite
    /// fibre, an unordered or non-positive stretch domain, or an invalid force-length family.
    pub fn new(
        t_max_pa: f64,
        force_length: ForceLength,
        fibre: [f64; 3],
        stretch_domain: (f64, f64),
    ) -> Result<Self, MaterialError> {
        let law = Self {
            t_max_pa,
            force_length,
            fibre,
            stretch_domain,
        };
        law.validate()?;
        Ok(law)
    }

    fn validate(&self) -> Result<(), MaterialError> {
        let bad = |what: &str| Err(MaterialError::Parameters { what: what.to_string() });
        if !(self.t_max_pa.is_finite() && self.t_max_pa > 0.0) {
            return bad("t_max_pa must be finite and positive");
        }
        let norm = fs_math::det::sqrt(self.fibre.iter().map(|x| x * x).sum());
        if !(self.fibre.iter().all(|x| x.is_finite()) && (norm - 1.0).abs() <= 1e-12) {
            return bad("the fibre must be a finite unit vector (|norm - 1| <= 1e-12)");
        }
        let (lo, hi) = self.stretch_domain;
        if !(lo.is_finite() && hi.is_finite() && 0.0 < lo && lo < hi) {
            return bad("the stretch domain must satisfy 0 < lo < hi, both finite");
        }
        self.force_length.validate()
    }
}

/// A passive hyperelastic law plus an active fibre stress along a fibre shared with it.
#[derive(Debug, Clone, PartialEq)]
pub struct ActiveTissue {
    /// The passive law.
    pub passive: Hyperelastic,
    /// The prescribed active stress.
    pub active: ActiveFibreStress,
}

fn det3<T: Real>(f: &[T; 9]) -> T {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

impl ActiveTissue {
    /// Compose a passive law and an active stress.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] when the active law is invalid, or when the passive law is an HGO card
    /// with fibre families and the active fibre is none of them (the fibre is shared with the passive card).
    pub fn new(passive: Hyperelastic, active: ActiveFibreStress) -> Result<Self, MaterialError> {
        active.validate()?;
        if let HyperelasticModel::HolzapfelGasserOgden { fibres, .. } = &passive.model {
            let shared = fibres.is_empty()
                || fibres
                    .iter()
                    .any(|a| (0..3).all(|k| (a[k] - active.fibre[k]).abs() <= 1e-12));
            if !shared {
                return Err(MaterialError::Parameters {
                    what: "the active fibre must be one of the passive HGO fibre families".to_string(),
                });
            }
        }
        Ok(Self { passive, active })
    }

    fn check_activation(activation: f64) -> Result<(), MaterialError> {
        if activation.is_finite() && (0.0..=1.0).contains(&activation) {
            Ok(())
        } else {
            Err(MaterialError::State {
                what: format!("activation {activation} outside [0, 1]"),
            })
        }
    }

    /// The fibre stretch `|F a0|` and its check against the declared domain (eight ULPs of inclusive
    /// roundoff, the same convention as the passive principal-stretch refusal).
    fn fibre_stretch(&self, f: &[f64; 9]) -> Result<f64, MaterialError> {
        let a = self.active.fibre;
        let fa: [f64; 3] =
            core::array::from_fn(|r| (0..3).fold(0.0, |s, c| s + f[3 * r + c] * a[c]));
        let lam = fs_math::det::sqrt(fa.iter().map(|x| x * x).sum());
        let (lo, hi) = self.active.stretch_domain;
        let roundoff = 8.0 * f64::EPSILON;
        if !lam.is_finite() || lam < lo * (1.0 - roundoff) || lam > hi * (1.0 + roundoff) {
            return Err(MaterialError::State {
                what: format!("fibre stretch {lam} outside the force-length domain [{lo}, {hi}]"),
            });
        }
        Ok(lam)
    }

    fn active_piola_generic<T: Real>(&self, f: &[T; 9], activation: f64) -> [T; 9] {
        let a0 = self.active.fibre;
        let fa: [T; 3] = core::array::from_fn(|r| {
            (0..3).fold(T::zero(), |s, c| s + f[3 * r + c] * T::from_f64(a0[c]))
        });
        let lam2 = fa[0] * fa[0] + fa[1] * fa[1] + fa[2] * fa[2];
        let lam = lam2.sqrt();
        let scale = T::from_f64(self.active.t_max_pa * activation)
            * self.active.force_length.factor(lam)
            * det3(f)
            * lam2.recip();
        core::array::from_fn(|ij| scale * fa[ij / 3] * T::from_f64(a0[ij % 3]))
    }

    /// The active first Piola-Kirchhoff stress alone (row-major 9-vector).
    ///
    /// # Errors
    /// [`MaterialError::State`] for an activation outside `[0, 1]`, a fibre stretch outside the domain, or an
    /// inadmissible deformation (`det F <= 0`).
    pub fn active_piola(&self, f: &[f64; 9], activation: f64) -> Result<[f64; 9], MaterialError> {
        Self::check_activation(activation)?;
        if !f.iter().all(|x| x.is_finite()) || !(det3(f) > 0.0) {
            return Err(MaterialError::State {
                what: format!("det F = {} (orientation lost or degenerate)", det3(f)),
            });
        }
        self.fibre_stretch(f)?;
        Ok(self.active_piola_generic(f, activation))
    }

    /// The total first Piola-Kirchhoff stress, passive plus active. At zero activation the passive law is
    /// returned unchanged, bit for bit.
    ///
    /// # Errors
    /// As [`Self::active_piola`], plus the passive law's own refusals.
    pub fn piola(&self, f: &[f64; 9], activation: f64) -> Result<[f64; 9], MaterialError> {
        Self::check_activation(activation)?;
        let mut p = self.passive.piola(f)?;
        if activation == 0.0 {
            return Ok(p);
        }
        let active = self.active_piola(f, activation)?;
        for (x, y) in p.iter_mut().zip(active) {
            *x += y;
        }
        Ok(p)
    }

    /// The consistent tangent `d P / d F` at fixed activation: the symmetric passive tangent plus the exact
    /// (forward-mode) derivative of the active stress, which is not symmetric.
    ///
    /// # Errors
    /// As [`Self::piola`].
    pub fn tangent(&self, f: &[f64; 9], activation: f64) -> Result<[[f64; 9]; 9], MaterialError> {
        Self::check_activation(activation)?;
        let mut a = self.passive.tangent(f)?;
        if activation == 0.0 {
            return Ok(a);
        }
        self.active_piola(f, activation)?;
        let vars: [Dual64<9>; 9] = core::array::from_fn(|i| Dual64::variable(f[i], i));
        let p = self.active_piola_generic(&vars, activation);
        for (row, p_row) in p.iter().enumerate() {
            for col in 0..9 {
                a[row][col] += p_row.eps[col];
            }
        }
        Ok(a)
    }

    /// The active power per reference volume, `P_a : Fdot` [W/m^3]: the player's mechanical input, booked
    /// separately from the passive stored energy.
    ///
    /// # Errors
    /// As [`Self::active_piola`].
    pub fn active_power(
        &self,
        f: &[f64; 9],
        f_dot: &[f64; 9],
        activation: f64,
    ) -> Result<f64, MaterialError> {
        let p = self.active_piola(f, activation)?;
        Ok(p.iter().zip(f_dot).map(|(a, b)| a * b).sum())
    }

    /// Thermodynamic declarations: the combination has no stored energy of its own (the passive part does),
    /// the active input is external work, and the active tangent is not symmetric.
    #[must_use]
    pub fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: false,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: false,
            failure_envelope: "none declared (prescribed active stress; refuses outside the force-length domain)",
        }
    }
}
