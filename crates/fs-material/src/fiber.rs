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
    /// Loading direction (+1/−1); zero until the first nonzero strain.
    pub dir: f64,
    /// Last committed strain.
    pub eps_prev: f64,
    /// Last committed stress.
    pub sig_prev: f64,
    /// Compression-side extremum recorded at reversals, seeded at −f_y/E₀.
    pub eps_min: f64,
    /// Tension-side extremum recorded at reversals, seeded at f_y/E₀.
    pub eps_max: f64,
    /// Maximum normalized excursion ξ seen so far (diagnostic only).
    /// The current branch curvature uses its own directional excursion.
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
        if step == 0.0 || step.signum().to_bits() == s.dir.to_bits() {
            return None;
        }
        let dir = step.signum();
        if s.dir == 0.0 {
            // Selecting the virgin loading direction is not a reversal: no
            // plastic excursion has occurred, so both signs retain R0 and
            // zero excursion history. Even a tiny nonzero first step selects
            // its branch; only an exact zero hold leaves it uninitialized.
            return Some(MpState {
                eps_0: dir * self.yield_strain(),
                sig_0: dir * self.fy,
                dir,
                ..s.clone()
            });
        }
        // New branch from the reversal point toward the opposite yield
        // asymptote: σ = b·E₀·ε + d·f_y(1−b). Intersect with the elastic
        // line through (ε_r, σ_r) of slope E₀.
        let (eps_r, sig_r) = (s.eps_prev, s.sig_prev);
        let eps_0 =
            (dir * self.fy * (1.0 - self.b) - sig_r + self.e0 * eps_r) / (self.e0 * (1.0 - self.b));
        let sig_0 = sig_r + self.e0 * (eps_0 - eps_r);
        let eps_min = s.eps_min.min(eps_r);
        let eps_max = s.eps_max.max(eps_r);
        // Menegotto–Pinto's plastic excursion is measured from the previous
        // extremum in the NEW loading direction to the NEW asymptote
        // intersection (the no-isotropic-shift Steel02 recurrence). Using
        // the old intersection instead spuriously softens subyield reversals;
        // reusing the largest excursion erases nested-loop directional memory.
        let eps_target = if dir > 0.0 { eps_max } else { eps_min };
        let xi = ((eps_target - eps_0) / self.yield_strain()).abs();
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
            eps_min,
            eps_max,
            xi_max: xi.max(s.xi_max),
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
        // The positive branch coordinates supply the initial tangent at
        // zero; the first nonzero update selects either virgin direction.
        MpState {
            eps_r: 0.0,
            sig_r: 0.0,
            eps_0: self.yield_strain(),
            sig_0: self.fy,
            r: self.r0,
            dir: 0.0,
            eps_prev: 0.0,
            sig_prev: 0.0,
            eps_min: -self.yield_strain(),
            eps_max: self.yield_strain(),
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
/// [`Self::response`] admits the parameters, committed history and the inclusive
/// ultimate-strain bound before returning a proposed state. At a closed tensile
/// gap it selects the compression-side tangent; this is a one-sided derivative,
/// not tensile stiffness.
///
/// Envelope: Mander, Priestley and Park (1988), Eqs. (3), (4), (6), (8),
/// doi:10.1061/(ASCE)0733-9445(1988)114:8(1804). The elastic cyclic lines here
/// remain a simplified rule, not that paper's full unloading/reloading model.
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

/// Admitted concrete response and its proposed, independently committable history.
#[derive(Debug, Clone, PartialEq)]
pub struct ManderResponse {
    /// Compression-positive stress in Pa; zero in an open tensile gap.
    pub stress_pa: f64,
    /// `d(stress)/d(compression strain)` in Pa. At gap closure the
    /// compression-side derivative is selected; post-peak tangents may be negative.
    pub tangent_pa: f64,
    /// History after accepting this strain. Evaluating a trial never mutates input.
    pub state: ManderState,
}

impl ManderConcrete {
    /// Construct with checks (requires Ec > Esec so the envelope shape
    /// parameter r is positive).
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on inadmissible parameters.
    pub fn new(fcc: f64, eps_cc: f64, ec: f64, eps_cu: f64) -> Result<Self, MaterialError> {
        if ![fcc, eps_cc, ec, eps_cu]
            .into_iter()
            .all(|value| value.is_finite() && value > 0.0)
            || eps_cu <= eps_cc
        {
            return Err(MaterialError::Parameters {
                what: "Mander needs finite positive fcc, eps_cc, Ec and eps_cu > eps_cc".to_string(),
            });
        }
        let esec = fcc / eps_cc;
        let r = ec / (ec - esec);
        if !(esec.is_finite() && esec > 0.0 && ec > esec && r.is_finite() && r > 1.0) {
            return Err(MaterialError::Parameters {
                what: format!(
                    "Mander Ec {ec} must exceed finite positive Esec {esec}, with representable r > 1"
                ),
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

    /// The analytic monotonic envelope and its derivative, compression positive.
    /// At zero strain the compression-side derivative is `Ec`; tensile strains
    /// have zero stress and stiffness. Use [`Self::response`] to enforce the
    /// original parameters, history and ultimate-strain domain.
    #[must_use]
    pub fn envelope(&self, strain: f64) -> (f64, f64) {
        if strain < 0.0 {
            return (0.0, 0.0); // no tension capacity
        }
        if strain == 0.0 {
            return (0.0, self.ec);
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

    /// Checked compression trial with elastic unloading and a retained tensile gap.
    ///
    /// The declared `eps_cu` is an inclusive validity bound, not a post-crushing
    /// damage law. A trial beyond it refuses; no stress clipping or history reset
    /// extends that domain. The public parameters and history are readmitted on
    /// every call, including during unloading. A retained maximum must lie on
    /// this same law's envelope within floating-point roundoff.
    pub fn response(
        &self,
        strain: f64,
        state: &ManderState,
    ) -> Result<ManderResponse, MaterialError> {
        Self::new(self.fcc, self.eps_cc, self.ec, self.eps_cu)?;
        let invalid = |what: &str| MaterialError::State {
            what: format!("Mander concrete: {what}"),
        };
        if !strain.is_finite() || strain > self.eps_cu {
            return Err(invalid("trial strain must be finite and at or below eps_cu"));
        }
        if !(state.eps_max.is_finite()
            && (0.0..=self.eps_cu).contains(&state.eps_max)
            && state.sig_max.is_finite()
            && state.sig_max >= 0.0)
        {
            return Err(invalid("committed maximum strain/stress is inadmissible"));
        }
        let expected_max = self.envelope(state.eps_max).0;
        let history_scale = expected_max.abs().max(state.sig_max);
        if !expected_max.is_finite()
            || (state.sig_max - expected_max).abs() > 64.0 * f64::EPSILON * history_scale
        {
            return Err(invalid(
                "committed maximum does not belong to the original envelope",
            ));
        }
        let (stress_pa, tangent_pa) = if strain >= state.eps_max {
            self.envelope(strain)
        } else {
            let eps_p = self.eps_plastic(state);
            if strain < eps_p {
                (0.0, 0.0)
            } else {
                let tangent = state.sig_max / (state.eps_max - eps_p);
                (tangent * (strain - eps_p), tangent)
            }
        };
        if !(stress_pa.is_finite() && stress_pa >= 0.0 && tangent_pa.is_finite()) {
            return Err(invalid("stress or tangent is unrepresentable"));
        }
        Ok(ManderResponse {
            stress_pa,
            tangent_pa,
            state: if strain > state.eps_max {
                ManderState {
                    eps_max: strain,
                    sig_max: stress_pa,
                }
            } else {
                state.clone()
            },
        })
    }
}

impl Uniaxial for ManderConcrete {
    type State = ManderState;

    fn initial_state(&self) -> ManderState {
        ManderState::default()
    }

    fn stress(&self, strain: f64, state: &ManderState) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |response| response.stress_pa)
    }

    fn tangent(&self, strain: f64, state: &ManderState) -> f64 {
        self.response(strain, state)
            .map_or(f64::NAN, |response| response.tangent_pa)
    }

    fn update_state(&self, strain: f64, state: &ManderState) -> ManderState {
        // This infallible trait cannot return a refusal. Preserve the accepted
        // history; callers that need the diagnosis use the checked trial owner.
        self.response(strain, state)
            .map_or_else(|_| state.clone(), |response| response.state)
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

// ------------------------------------------------------------ Wool felt

/// Wool-felt compressive hysteresis (music bead
/// `frankensim-music-t-piano-felt-87zbd`): the piano-hammer felt law as a
/// NEW implementor of the existing [`Uniaxial`] contract — never Mander
/// concrete or Menegotto–Pinto steel wearing a felt hat (their envelopes
/// soften; felt STIFFENS, and its residual state is compaction, not a
/// modulus-derived plastic offset).
///
/// Physics, rate-independent half only:
/// - Virgin loading follows the stiffening power law the felt literature
///   established for hammer compression, `σ = f_ref · (ε/ε_ref)^p` with
///   `p > 1` (measured piano-hammer exponents cluster ≈ 2.2–3.5).
/// - Unloading/reloading inside a past excursion runs on a steeper
///   power-law path anchored between the crush point `(ε_r, 0)` and the
///   excursion peak `(ε_max, σ_max)`: `σ = σ_max · x^q`,
///   `x = (ε−ε_r)/(ε_max−ε_r)`, with `q ≥ p`.
/// - Residual crush is a fraction of the excursion: `ε_r = c · ε_max`
///   (felt compacts; hammers need voicing).
/// - Zero tension capacity.
///
/// Dissipation is guaranteed BY CONSTRUCTION: over a virgin excursion to
/// `ε_max` and return, the envelope work is
/// `W_env = σ_max · ε_max / (p+1)` and the unload path returns
/// `W_ret = σ_max · (ε_max − ε_r) / (q+1)`, so
/// `W_env / W_ret = (q+1) / ((p+1)(1−c)) > 1` whenever `q ≥ p` and
/// `c > 0` — both enforced at construction, so the mt-004 loop-area gate
/// holds structurally, not accidentally.
///
/// The RATE-DEPENDENT half of felt (loading-rate stiffening, sub-loop
/// creep hysteresis) is deliberately NOT here: it belongs to the
/// [`crate::visco::FractionalZener`] → `lower_to_prony` island the crate
/// already names for felt. Closed sub-loops inside a committed excursion
/// are elastic in THIS law; their loss is the Prony island's claim.
#[derive(Debug, Clone, PartialEq)]
pub struct WoolFelt {
    /// Stress at the reference strain (Pa, compression positive).
    pub f_ref: f64,
    /// Reference compressive strain (dimensionless engineering strain).
    pub eps_ref: f64,
    /// Loading (envelope) exponent `p > 1`.
    pub p: f64,
    /// Unload/reload exponent `q ≥ p`.
    pub q: f64,
    /// Residual-crush fraction `c ∈ (0, 1)` of the maximum excursion.
    /// `c = 0` with `q = p` is the conservative power-law limit.
    pub crush_fraction: f64,
    /// Densification strain — validity bound where felt bottoms out.
    pub eps_densify: f64,
}

/// Wool-felt committed state.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct WoolFeltState {
    /// Maximum compressive strain ever committed.
    pub eps_max: f64,
    /// Envelope stress at that excursion.
    pub sig_max: f64,
}

impl WoolFelt {
    /// Construct with admission checks.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] unless `f_ref > 0`, `eps_ref > 0`,
    /// `p > 1`, `q ≥ p`, `0 < crush_fraction < 1` (or `c = 0, q = p`), and
    /// `eps_densify > eps_ref` — the parameter region where the
    /// dissipation proof above holds strictly.
    pub fn new(
        f_ref: f64,
        eps_ref: f64,
        p: f64,
        q: f64,
        crush_fraction: f64,
        eps_densify: f64,
    ) -> Result<Self, MaterialError> {
        if !(f_ref > 0.0
            && f_ref.is_finite()
            && eps_ref > 0.0
            && p > 1.0
            && q >= p
            && q.is_finite()
            && (crush_fraction > 0.0
                || (crush_fraction == 0.0
                    && q.partial_cmp(&p) == Some(core::cmp::Ordering::Equal)))
            && crush_fraction < 1.0
            && eps_densify > eps_ref)
        {
            return Err(MaterialError::Parameters {
                what: "wool felt needs f_ref>0, eps_ref>0, p>1, q>=p, 0<crush<1 or crush=0 with q=p, densify>ref"
                    .to_string(),
            });
        }
        Ok(Self {
            f_ref,
            eps_ref,
            p,
            q,
            crush_fraction,
            eps_densify,
        })
    }

    /// Virgin envelope `(σ, dσ/dε)` at compressive strain `ε ≥ 0`.
    #[must_use]
    pub fn envelope(&self, strain: f64) -> (f64, f64) {
        if strain <= 0.0 {
            return (0.0, 0.0); // zero tension capacity
        }
        let x = strain / self.eps_ref;
        let sig = self.f_ref * det::pow(x, self.p);
        let dsig = self.f_ref * self.p * det::pow(x, self.p - 1.0) / self.eps_ref;
        (sig, dsig)
    }

    /// Residual (crush) strain for a committed state.
    #[must_use]
    pub fn eps_residual(&self, state: &WoolFeltState) -> f64 {
        self.crush_fraction * state.eps_max
    }
}

impl Uniaxial for WoolFelt {
    type State = WoolFeltState;

    fn initial_state(&self) -> WoolFeltState {
        WoolFeltState::default()
    }

    fn stress(&self, strain: f64, state: &WoolFeltState) -> f64 {
        if strain >= state.eps_max {
            return self.envelope(strain).0;
        }
        let eps_r = self.eps_residual(state);
        if strain <= eps_r {
            return 0.0;
        }
        let x = (strain - eps_r) / (state.eps_max - eps_r);
        state.sig_max * det::pow(x, self.q)
    }

    fn tangent(&self, strain: f64, state: &WoolFeltState) -> f64 {
        if strain >= state.eps_max {
            return self.envelope(strain).1;
        }
        let eps_r = self.eps_residual(state);
        if strain <= eps_r {
            return 0.0;
        }
        let span = state.eps_max - eps_r;
        let x = (strain - eps_r) / span;
        state.sig_max * self.q * det::pow(x, self.q - 1.0) / span
    }

    fn update_state(&self, strain: f64, state: &WoolFeltState) -> WoolFeltState {
        if strain > state.eps_max {
            let (sig, _) = self.envelope(strain);
            WoolFeltState {
                eps_max: strain,
                sig_max: sig,
            }
        } else {
            state.clone()
        }
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true, // (q+1) >= (p+1)(1-c) by construction
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "densification beyond eps_densify (validity bound); zero tension",
        }
    }

    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.wool-felt",
            "0.1.0",
            Ambition::Solid,
            vec![
                "uniaxial fiber, compression positive".to_string(),
                "stiffening power-law envelope (hammer-felt literature form)".to_string(),
                "power-law unload through fractional residual crush".to_string(),
                "rate independence; loading-rate and sub-loop loss belong to the \
                 FractionalZener/Prony island"
                    .to_string(),
            ],
            ValidityDomain::unconstrained().with("strain", 0.0, self.eps_densify),
            vec![
                "zero tension capacity".to_string(),
                "no rate dependence here (visco island's claim)".to_string(),
                "no densification hardening beyond eps_densify".to_string(),
                "closed sub-loops inside a committed excursion are elastic".to_string(),
            ],
            0.15,
        )
    }
}
