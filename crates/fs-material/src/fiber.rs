//! Uniaxial fiber laws — the reinforced-concrete flagship pair: Mander
//! confined concrete + Menegotto–Pinto steel. These feed fiber-section
//! beams; their tangent contract is identical to the 3D laws (consistent
//! dσ/dε at the same committed-state update, FD-gated).

use crate::{MaterialAdmissibility, MaterialError};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

/// A uniaxial (fiber) constitutive law with hysteretic state.
pub trait Uniaxial {
    /// Committed internal variables.
    type State: Clone + PartialEq + core::fmt::Debug;

    /// The virgin state.
    fn initial_state(&self) -> Self::State;

    /// Stress at `strain` given the last committed state.
    fn stress(&self, strain: f64, state: &Self::State) -> f64;

    /// Consistent tangent dσ/dε at the same update.
    fn tangent(&self, strain: f64, state: &Self::State) -> f64;

    /// Commit the state at `strain`.
    fn update_state(&self, strain: f64, state: &Self::State) -> Self::State;

    /// Thermodynamic declarations.
    fn admissibility(&self) -> MaterialAdmissibility;

    /// The model card.
    fn card(&self) -> ModelCard;
}

// ------------------------------------------------------- Menegotto–Pinto

/// Menegotto–Pinto steel with Bauschinger-effect curvature degradation
/// (R0, a1, a2 per the classic formulation).
#[derive(Debug, Clone, PartialEq)]
pub struct MenegottoPintoSteel {
    /// Elastic modulus E₀ (Pa).
    pub e0: f64,
    /// Yield stress f_y (Pa).
    pub fy: f64,
    /// Hardening ratio b (post-yield slope = b·E₀).
    pub b: f64,
    /// Initial curvature parameter R₀ (≈ 20).
    pub r0: f64,
    /// Curvature degradation a₁ (≈ 18.5).
    pub a1: f64,
    /// Curvature degradation a₂ (≈ 0.15).
    pub a2: f64,
}

/// Menegotto–Pinto branch state.
#[derive(Debug, Clone, PartialEq)]
pub struct MpState {
    /// Strain at the last reversal.
    pub eps_r: f64,
    /// Stress at the last reversal.
    pub sig_r: f64,
    /// Asymptote-intersection strain of the current branch.
    pub eps_0: f64,
    /// Asymptote-intersection stress of the current branch.
    pub sig_0: f64,
    /// Current branch curvature R.
    pub r: f64,
    /// Loading direction of the current branch (+1/−1).
    pub dir: f64,
    /// Last committed strain.
    pub eps_prev: f64,
    /// Last committed stress.
    pub sig_prev: f64,
    /// Maximum normalized excursion ξ (drives R degradation).
    pub xi_max: f64,
}

impl MenegottoPintoSteel {
    /// Construct with checks.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on non-positive E₀/f_y or b ∉ [0, 1).
    pub fn new(e0: f64, fy: f64, b: f64) -> Result<Self, MaterialError> {
        if !(e0 > 0.0 && fy > 0.0 && (0.0..1.0).contains(&b)) {
            return Err(MaterialError::Parameters {
                what: format!("E0 {e0}, fy {fy}, b {b} out of range"),
            });
        }
        Ok(MenegottoPintoSteel {
            e0,
            fy,
            b,
            r0: 20.0,
            a1: 18.5,
            a2: 0.15,
        })
    }

    fn yield_strain(&self) -> f64 {
        self.fy / self.e0
    }

    /// Normalized branch coordinates and curve evaluation.
    fn curve(&self, strain: f64, s: &MpState) -> (f64, f64) {
        let de = s.eps_0 - s.eps_r;
        let ds = s.sig_0 - s.sig_r;
        let eps_star = (strain - s.eps_r) / de;
        let abs_pow = det::pow(eps_star.abs().max(1e-300), s.r);
        let denom = det::pow(1.0 + abs_pow, 1.0 / s.r);
        let sig_star = self.b * eps_star + (1.0 - self.b) * eps_star / denom;
        // dσ*/dε* = b + (1−b)(1+|ε*|^R)^(−1/R−1)
        let dsig_star = self.b + (1.0 - self.b) * det::pow(1.0 + abs_pow, -1.0 / s.r - 1.0);
        (s.sig_r + sig_star * ds, dsig_star * ds / de)
    }

    fn reversed(&self, strain: f64, s: &MpState) -> Option<MpState> {
        let step = strain - s.eps_prev;
        if step.abs() < f64::MIN_POSITIVE || step.signum().to_bits() == s.dir.to_bits() {
            return None;
        }
        // New branch from the reversal point toward the opposite yield
        // asymptote: σ = b·E₀·ε + d·f_y(1−b). Intersect with the elastic
        // line through (ε_r, σ_r) of slope E₀.
        let dir = step.signum();
        let (eps_r, sig_r) = (s.eps_prev, s.sig_prev);
        let eps_0 =
            (dir * self.fy * (1.0 - self.b) - sig_r + self.e0 * eps_r) / (self.e0 * (1.0 - self.b));
        let sig_0 = sig_r + self.e0 * (eps_0 - eps_r);
        let xi = ((eps_r - s.eps_0) / self.yield_strain())
            .abs()
            .max(s.xi_max);
        let r = self.r0 - self.a1 * xi / (self.a2 + xi);
        Some(MpState {
            eps_r,
            sig_r,
            eps_0,
            sig_0,
            r,
            dir,
            eps_prev: s.eps_prev,
            sig_prev: s.sig_prev,
            xi_max: xi,
        })
    }

    fn effective(&self, strain: f64, state: &MpState) -> MpState {
        self.reversed(strain, state)
            .unwrap_or_else(|| state.clone())
    }
}

impl Uniaxial for MenegottoPintoSteel {
    type State = MpState;

    fn initial_state(&self) -> MpState {
        // Virgin monotonic branch toward +yield (mirrored on first
        // negative excursion by the reversal logic with eps_prev = 0).
        MpState {
            eps_r: 0.0,
            sig_r: 0.0,
            eps_0: self.yield_strain(),
            sig_0: self.fy,
            r: self.r0,
            dir: 1.0,
            eps_prev: 0.0,
            sig_prev: 0.0,
            xi_max: 0.0,
        }
    }

    fn stress(&self, strain: f64, state: &MpState) -> f64 {
        let s = self.effective(strain, state);
        self.curve(strain, &s).0
    }

    fn tangent(&self, strain: f64, state: &MpState) -> f64 {
        let s = self.effective(strain, state);
        self.curve(strain, &s).1
    }

    fn update_state(&self, strain: f64, state: &MpState) -> MpState {
        let mut s = self.effective(strain, state);
        let (sig, _) = self.curve(strain, &s);
        s.eps_prev = strain;
        s.sig_prev = sig;
        s
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: false,
            dissipation_nonnegative: true, // closed loops dissipate
            polyconvex: None,
            tangent_symmetric: true, // scalar
            failure_envelope: "no rupture criterion (strain validity bound on the card)",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.menegotto-pinto-steel",
            "0.1.0",
            Ambition::Solid,
            vec![
                "uniaxial fiber".to_string(),
                "rate independence".to_string(),
                "Bauschinger via R degradation (R0/a1/a2)".to_string(),
            ],
            ValidityDomain::unconstrained().with("strain-magnitude", 0.0, 0.08),
            vec![
                "no buckling of the bar".to_string(),
                "no low-cycle fatigue rupture".to_string(),
            ],
            0.05,
        )
    }
}

// ------------------------------------------------------------- Mander

/// Mander confined-concrete envelope with elastic unload/reload lines
/// (compression positive; zero tension capacity — documented no-claims).
#[derive(Debug, Clone, PartialEq)]
pub struct ManderConcrete {
    /// Confined compressive strength f′cc (Pa, positive).
    pub fcc: f64,
    /// Strain at f′cc.
    pub eps_cc: f64,
    /// Initial (elastic) modulus Ec (Pa).
    pub ec: f64,
    /// Ultimate (crushing) strain — validity bound.
    pub eps_cu: f64,
}

/// Mander committed state.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ManderState {
    /// Maximum compressive strain ever committed.
    pub eps_max: f64,
    /// Stress at that excursion (envelope value).
    pub sig_max: f64,
}

impl ManderConcrete {
    /// Construct with checks (requires Ec > Esec so the envelope shape
    /// parameter r is positive).
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on inadmissible parameters.
    pub fn new(fcc: f64, eps_cc: f64, ec: f64, eps_cu: f64) -> Result<Self, MaterialError> {
        if !(fcc > 0.0 && eps_cc > 0.0 && ec > 0.0 && eps_cu > eps_cc) {
            return Err(MaterialError::Parameters {
                what: "need fcc, eps_cc, Ec > 0 and eps_cu > eps_cc".to_string(),
            });
        }
        let esec = fcc / eps_cc;
        if ec <= esec {
            return Err(MaterialError::Parameters {
                what: format!("Ec {ec} must exceed the secant modulus {esec}"),
            });
        }
        Ok(ManderConcrete {
            fcc,
            eps_cc,
            ec,
            eps_cu,
        })
    }

    fn shape_r(&self) -> f64 {
        let esec = self.fcc / self.eps_cc;
        self.ec / (self.ec - esec)
    }

    /// The monotonic envelope σ(ε) and its derivative (ε ≥ 0,
    /// compression positive).
    #[must_use]
    pub fn envelope(&self, strain: f64) -> (f64, f64) {
        if strain <= 0.0 {
            return (0.0, 0.0); // no tension capacity
        }
        let r = self.shape_r();
        let x = strain / self.eps_cc;
        let xr = det::pow(x, r);
        let denom = r - 1.0 + xr;
        let sig = self.fcc * x * r / denom;
        let dsig_dx = self.fcc * r * (r - 1.0) * (1.0 - xr) / (denom * denom);
        (sig, dsig_dx / self.eps_cc)
    }

    /// Plastic (residual) strain after unloading from `state.eps_max`.
    fn eps_plastic(&self, state: &ManderState) -> f64 {
        (state.eps_max - state.sig_max / self.ec).max(0.0)
    }
}

impl Uniaxial for ManderConcrete {
    type State = ManderState;

    fn initial_state(&self) -> ManderState {
        ManderState::default()
    }

    fn stress(&self, strain: f64, state: &ManderState) -> f64 {
        if strain >= state.eps_max {
            return self.envelope(strain).0;
        }
        // Unload/reload line between (eps_p, 0) and (eps_max, sig_max).
        let eps_p = self.eps_plastic(state);
        if strain <= eps_p {
            return 0.0;
        }
        state.sig_max * (strain - eps_p) / (state.eps_max - eps_p)
    }

    fn tangent(&self, strain: f64, state: &ManderState) -> f64 {
        if strain >= state.eps_max {
            return self.envelope(strain).1;
        }
        let eps_p = self.eps_plastic(state);
        if strain <= eps_p {
            return 0.0;
        }
        state.sig_max / (state.eps_max - eps_p)
    }

    fn update_state(&self, strain: f64, state: &ManderState) -> ManderState {
        if strain > state.eps_max {
            let (sig, _) = self.envelope(strain);
            ManderState {
                eps_max: strain,
                sig_max: sig,
            }
        } else {
            state.clone()
        }
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: false,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "crushing at eps_cu (validity bound); zero tension",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.mander-concrete",
            "0.1.0",
            Ambition::Solid,
            vec![
                "uniaxial fiber, compression positive".to_string(),
                "confinement folded into (fcc, eps_cc)".to_string(),
                "elastic unload/reload lines (simplified cyclic rule)".to_string(),
            ],
            ValidityDomain::unconstrained().with("strain", 0.0, self.eps_cu),
            vec![
                "zero tension capacity (no tension stiffening)".to_string(),
                "no cyclic stiffness degradation beyond the residual-strain rule".to_string(),
            ],
            0.10,
        )
    }
}
