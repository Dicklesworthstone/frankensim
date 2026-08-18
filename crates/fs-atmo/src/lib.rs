//! fs-atmo — atmosphere foundation (L2). Bead frankensim-wf-root-guzez.4.5
//! (E3.3a, Wright Flyer program).
//!
//! Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md
//! §5.4 (ROUND 6 steady state); evidence: air-state-v1.json (E1.8) under the
//! frozen registry (E1.7). Conventions: frd/NED (frame-conventions-v1) —
//! +z DOWN; heights here are ALTITUDE h above the aerodynamic ground plane
//! (h = −z), converted at the API boundary.
//!
//! E3.3a scope (V-04a, analytic construction):
//! - `FlatSiteLogLaw` — ONE scenario-level effective z₀ (pointwise z₀
//!   insertion is FORBIDDEN while solenoidal claims stand: a fetch-varying
//!   z₀ in U(x, h; z₀(x))·e_x makes ∂U/∂x ≠ 0 — Round-2 fix).
//! - Wall-compatible vector-potential turbulence: u = ∇×ψ with the
//!   horizontal potential components carrying sin(k_h·h) vertical shape, so
//!   u_vertical(h = 0) = 0 IDENTICALLY (bitwise) and ∇·u ≡ 0 analytically.
//! - Exact analytic derivatives (the gradient is derived term-wise, never
//!   finite-differenced internally).
//! - `sample_air_state` — velocity + gradient + ρ/μ/T/p with provenance;
//!   Re and q derive from the SAME provenance-bound state.
//! - Deterministic philox stream partitioning: mode i draws from
//!   StreamKey{seed, ATMO_KERNEL, tile = i} — counter-addressed, so mode
//!   i's amplitudes are independent of how many other modes exist.
//!
//! Boundary: the exact-discrete OU amplitude EVOLUTION (sequential,
//! checkpointed state) and the Mann-class spectral-tensor FIT are E3.3b.
//! Here amplitudes are a static per-realization draw with the plan's
//! deterministic mean-advection phase φ_k = k·(x − U_adv·t); the declared
//! amplitude decay is an analytic-construction placeholder whose statistics
//! carry NO claim (V-04b territory).

use fs_math::det;
use fs_rand::StreamKey;

/// Registered fs-rand kernel id for atmosphere draws ("ATMO").
pub const ATMO_KERNEL: u32 = 0x41544D4F;
/// von Kármán constant (declared model constant).
pub const KAPPA: f64 = 0.40;
/// Admitted z₀ domain [m] (air-state-v1 prior support, widened one decade).
pub const MIN_Z0_M: f64 = 1.0e-5;
/// Upper z₀ bound [m].
pub const MAX_Z0_M: f64 = 1.0;
/// Mode-count cap (refusals at cap AND cap+1).
pub const MAX_MODES: usize = 256;
/// The fixed tick rate the phase clock uses [Hz].
pub const TICK_HZ: f64 = 120.0;

/// A typed refusal (workspace law).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// Human-readable diagnosis.
    pub message: String,
    /// Ranked repairs, most likely fix first.
    pub ranked_repairs: Vec<String>,
}

fn refuse(code: &'static str, message: String, repair: &str) -> Refusal {
    Refusal {
        code,
        message,
        ranked_repairs: vec![repair.into()],
    }
}

/// The historical 1903 mean-wind mode: one scenario-level effective z₀
/// over the certified launch region (plan §5.4).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlatSiteLogLaw {
    /// Scenario-level effective roughness length [m].
    pub scenario_effective_z0_m: f64,
    /// Displacement height d [m].
    pub displacement_height_m: f64,
    /// Reference height above ground [m] (where the reference speed holds).
    pub reference_height_m: f64,
    /// Mean speed at the reference height [m/s].
    pub reference_speed_mps: f64,
}

impl FlatSiteLogLaw {
    /// Validate the law's parameters.
    ///
    /// # Errors
    /// `non-finite-input`, `z0-outside-domain`, `displacement-invalid`,
    /// `reference-height-invalid`, `reference-speed-invalid`.
    pub fn admit(&self) -> Result<(), Refusal> {
        let all = [
            self.scenario_effective_z0_m,
            self.displacement_height_m,
            self.reference_height_m,
            self.reference_speed_mps,
        ];
        if !all.iter().all(|v| v.is_finite()) {
            return Err(refuse(
                "non-finite-input",
                format!("{self:?}"),
                "finite parameters only",
            ));
        }
        if !(MIN_Z0_M..=MAX_Z0_M).contains(&self.scenario_effective_z0_m) {
            return Err(refuse(
                "z0-outside-domain",
                format!(
                    "z0 {} m outside [{MIN_Z0_M}, {MAX_Z0_M}]",
                    self.scenario_effective_z0_m
                ),
                "draw z0 from the air-state-v1 log-uniform prior [3e-4, 1e-2]",
            ));
        }
        if self.displacement_height_m < 0.0 || self.displacement_height_m > 1.0 {
            return Err(refuse(
                "displacement-invalid",
                format!("d {} m outside [0, 1]", self.displacement_height_m),
                "the 1903 sand plain prior is [0, 0.05] m",
            ));
        }
        if self.reference_height_m <= self.displacement_height_m + self.scenario_effective_z0_m {
            return Err(refuse(
                "reference-height-invalid",
                format!(
                    "reference height {} m must exceed d + z0 = {} m",
                    self.reference_height_m,
                    self.displacement_height_m + self.scenario_effective_z0_m
                ),
                "the instrument-height priors (1.5-10 m) all satisfy this",
            ));
        }
        if self.reference_speed_mps < 0.0 {
            return Err(refuse(
                "reference-speed-invalid",
                format!("reference speed {} m/s negative", self.reference_speed_mps),
                "speeds are magnitudes; direction is the flow frame's job",
            ));
        }
        Ok(())
    }

    /// Friction velocity u* implied by the reference point [m/s].
    #[must_use]
    pub fn u_star(&self) -> f64 {
        let arg =
            (self.reference_height_m - self.displacement_height_m) / self.scenario_effective_z0_m;
        self.reference_speed_mps * KAPPA / det::ln(arg)
    }

    /// Mean speed at altitude h above the ground plane [m/s]:
    /// U(h) = (u*/κ)·ln((h − d)/z₀) for h − d ≥ z₀, else 0 (inside the
    /// roughness sublayer the log law does not apply; the mean is clamped
    /// to zero rather than extrapolated negative).
    #[must_use]
    pub fn speed(&self, h_m: f64) -> f64 {
        let excess = h_m - self.displacement_height_m;
        if excess <= self.scenario_effective_z0_m {
            return 0.0;
        }
        self.u_star() / KAPPA * det::ln(excess / self.scenario_effective_z0_m)
    }

    /// Analytic dU/dh [1/s]: u*/(κ·(h − d)) above the sublayer, else 0.
    #[must_use]
    pub fn dspeed_dh(&self, h_m: f64) -> f64 {
        let excess = h_m - self.displacement_height_m;
        if excess <= self.scenario_effective_z0_m {
            return 0.0;
        }
        self.u_star() / (KAPPA * excess)
    }
}

/// One turbulence mode of the wall-compatible vector potential.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Mode {
    /// Horizontal wavevector components [rad/m].
    kx: f64,
    ky: f64,
    /// Vertical wavenumber [rad/m].
    kh: f64,
    /// Potential amplitudes (a_x, a_y, a_z) [m²/s].
    a: [f64; 3],
    /// Phase offset [rad].
    phi: f64,
}

/// The E3.3a turbulence field: a static per-realization modal draw with
/// deterministic mean-advection phase (OU evolution is E3.3b).
#[derive(Clone, Debug, PartialEq)]
pub struct TurbulenceField {
    modes: Vec<Mode>,
    /// Advection speed for the frozen-phase clock [m/s].
    u_adv_mps: f64,
}

/// Sampled velocity + exact analytic gradient at one point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlowSample {
    /// Velocity (u_east-ish frame: [along-wind, cross-wind, VERTICAL-UP])
    /// — converted to frd/NED by the caller-facing air-state API.
    pub u: [f64; 3],
    /// Gradient du_i/dx_j, row i component, column j ∈ {x, y, h}.
    pub grad: [[f64; 3]; 3],
}

impl TurbulenceField {
    /// Build a field: `n_modes` philox-partitioned modes at `seed`, with
    /// declared analytic-construction amplitudes (rms `sigma_mps`, integral
    /// scale `length_m`) and the given advection speed.
    ///
    /// # Errors
    /// `mode-count-invalid` (0 or above [`MAX_MODES`], tested at cap AND
    /// cap+1); `non-finite-input`; `turbulence-params-invalid`.
    pub fn build(
        seed: u64,
        n_modes: usize,
        sigma_mps: f64,
        length_m: f64,
        u_adv_mps: f64,
    ) -> Result<TurbulenceField, Refusal> {
        if n_modes == 0 || n_modes > MAX_MODES {
            return Err(refuse(
                "mode-count-invalid",
                format!("{n_modes} modes outside [1, {MAX_MODES}]"),
                "64 modes is the browser-tier default",
            ));
        }
        if !(sigma_mps.is_finite() && length_m.is_finite() && u_adv_mps.is_finite()) {
            return Err(refuse(
                "non-finite-input",
                "sigma/length/u_adv".into(),
                "finite only",
            ));
        }
        if sigma_mps < 0.0 || length_m <= 0.0 {
            return Err(refuse(
                "turbulence-params-invalid",
                format!("sigma {sigma_mps}, length {length_m}"),
                "sigma >= 0, length > 0",
            ));
        }
        let mut modes = Vec::with_capacity(n_modes);
        let base_k = core::f64::consts::TAU / (8.0 * length_m);
        for i in 0..n_modes {
            // Counter-addressed partition: mode i's stream is a pure
            // function of (seed, ATMO_KERNEL, i) — independent of n_modes.
            let mut s = StreamKey {
                seed,
                kernel: ATMO_KERNEL,
                tile: i as u32,
            }
            .stream();
            // Log-spaced shell + random direction (deterministic draws).
            let shell = base_k * det::exp(det::ln(64.0) * s.next_f64());
            let angle = core::f64::consts::TAU * s.next_f64();
            let kx = shell * det::cos(angle);
            let ky = shell * det::sin(angle);
            let kh = shell * (0.5 + s.next_f64());
            // Declared decay: potential amplitude ∝ shell^(-11/6)·σ/√N —
            // an analytic-construction placeholder (statistics are E3.3b).
            let scale = sigma_mps / det::sqrt(n_modes as f64) / shell
                * det::exp(-11.0 / 6.0 * det::ln(shell / base_k));
            let a = [
                scale * s.next_normal(),
                scale * s.next_normal(),
                scale * s.next_normal(),
            ];
            let phi = core::f64::consts::TAU * s.next_f64();
            modes.push(Mode { kx, ky, kh, a, phi });
        }
        Ok(TurbulenceField { modes, u_adv_mps })
    }

    /// Sample velocity + exact analytic gradient at position (x, y, h) and
    /// `tick` (the frozen-phase clock t = tick/120 s).
    ///
    /// Construction (per mode, θ = kx·(x − U_adv·t) + ky·y + φ):
    ///   ψ_x = a_x·sin(k_h h)·cos θ, ψ_y = a_y·sin(k_h h)·cos θ,
    ///   ψ_z = a_z·cos(k_h h)·cos θ,  u = ∇×ψ.
    /// sin(k_h·0) = 0 makes the vertical component vanish at the wall
    /// IDENTICALLY, and a curl is divergence-free analytically.
    #[must_use]
    pub fn sample(&self, x_m: f64, y_m: f64, h_m: f64, tick: u64) -> FlowSample {
        let t = tick as f64 / TICK_HZ;
        let mut u = [0.0f64; 3];
        let mut grad = [[0.0f64; 3]; 3];
        for m in &self.modes {
            let theta = m.kx * (x_m - self.u_adv_mps * t) + m.ky * y_m + m.phi;
            let (st, ct) = (det::sin(theta), det::cos(theta));
            let (sh, ch) = (det::sin(m.kh * h_m), det::cos(m.kh * h_m));
            let [ax, ay, az] = m.a;
            // u = ∇×ψ (derived term-wise; see the battery's FD cross-check).
            let ux = -az * m.ky * ch * st - ay * m.kh * ch * ct;
            let uy = ax * m.kh * ch * ct + az * m.kx * ch * st;
            let uz = (ax * m.ky - ay * m.kx) * sh * st;
            u[0] += ux;
            u[1] += uy;
            u[2] += uz;
            // Analytic partials. d(st)/dx = kx·ct, d(ct)/dx = −kx·st, etc.
            grad[0][0] += -az * m.ky * m.kx * ch * ct + ay * m.kh * m.kx * ch * st;
            grad[0][1] += -az * m.ky * m.ky * ch * ct + ay * m.kh * m.ky * ch * st;
            grad[0][2] += az * m.ky * m.kh * sh * st + ay * m.kh * m.kh * sh * ct;
            grad[1][0] += ax * m.kh * (-m.kx) * ch * st + az * m.kx * m.kx * ch * ct;
            grad[1][1] += ax * m.kh * (-m.ky) * ch * st + az * m.kx * m.ky * ch * ct;
            grad[1][2] += -ax * m.kh * m.kh * sh * ct - az * m.kx * m.kh * sh * st;
            let c = ax * m.ky - ay * m.kx;
            grad[2][0] += c * m.kx * sh * ct;
            grad[2][1] += c * m.ky * sh * ct;
            grad[2][2] += c * m.kh * ch * st;
        }
        FlowSample { u, grad }
    }
}

/// Scenario air constants with provenance (E1.8 dossier values enter here).
#[derive(Clone, Debug, PartialEq)]
pub struct AirScenario {
    /// Density [kg/m³].
    pub rho_kg_m3: f64,
    /// Dynamic viscosity [kg/(m·s)].
    pub mu_kg_m_s: f64,
    /// Temperature [K].
    pub temperature_k: f64,
    /// Pressure [Pa].
    pub pressure_pa: f64,
    /// Provenance string (dossier record + derivation).
    pub provenance: &'static str,
}

/// The Dec-17 1903 ensemble-mean air state (air-state-v1 derivations).
pub const DEC17_AIR: AirScenario = AirScenario {
    rho_kg_m3: 1.294,
    mu_kg_m_s: 1.721e-5,
    temperature_k: 274.3,
    pressure_pa: 101_900.0,
    provenance: "air-state-v1 (E1.8): LSS log 34F/30.1inHg; rho=p/(RT); Sutherland mu",
};

/// One provenance-bound air-state sample (plan §5.4 API).
#[derive(Clone, Debug, PartialEq)]
pub struct AirState {
    /// Total velocity [m/s] in flow frame [along-wind, cross, vertical-up].
    pub velocity_mps: [f64; 3],
    /// Exact analytic gradient of the total velocity.
    pub grad: [[f64; 3]; 3],
    /// Density [kg/m³] (same provenance as the velocity's scenario).
    pub rho_kg_m3: f64,
    /// Dynamic viscosity [kg/(m·s)].
    pub mu_kg_m_s: f64,
    /// Temperature [K].
    pub temperature_k: f64,
    /// Pressure [Pa].
    pub pressure_pa: f64,
    /// Provenance of the air constants.
    pub provenance: &'static str,
}

impl AirState {
    /// Dynamic pressure q = ½ρ|V|² [Pa] — derives from the SAME state.
    #[must_use]
    pub fn dynamic_pressure_pa(&self) -> f64 {
        let v2: f64 = self.velocity_mps.iter().map(|v| v * v).sum();
        0.5 * self.rho_kg_m3 * v2
    }

    /// Reynolds number over `chord_m` — same-state ρ, μ, |V|.
    #[must_use]
    pub fn reynolds(&self, chord_m: f64) -> f64 {
        let v: f64 = det::sqrt(self.velocity_mps.iter().map(|v| v * v).sum());
        self.rho_kg_m3 * v * chord_m / self.mu_kg_m_s
    }
}

/// The E3.3a atmosphere: mean law + turbulence + air scenario.
#[derive(Clone, Debug, PartialEq)]
pub struct Atmosphere {
    /// Mean-wind law.
    pub mean: FlatSiteLogLaw,
    /// Turbulence field.
    pub turbulence: TurbulenceField,
    /// Air constants with provenance.
    pub air: AirScenario,
}

impl Atmosphere {
    /// Sample the full air state at (x, y, h) and tick.
    ///
    /// # Errors
    /// `non-finite-input`; `below-surface-query` (h < 0 refuses — the
    /// ground plane bounds the domain; contact is fs-flyer's job).
    pub fn sample_air_state(
        &self,
        x_m: f64,
        y_m: f64,
        h_m: f64,
        tick: u64,
    ) -> Result<AirState, Refusal> {
        if !(x_m.is_finite() && y_m.is_finite() && h_m.is_finite()) {
            return Err(refuse(
                "non-finite-input",
                format!("({x_m}, {y_m}, {h_m})"),
                "finite",
            ));
        }
        if h_m < 0.0 {
            return Err(refuse(
                "below-surface-query",
                format!("h = {h_m} m is below the aerodynamic ground plane"),
                "the FlatnessCertificate plane bounds the domain at h = 0",
            ));
        }
        self.mean.admit()?;
        let turb = self.turbulence.sample(x_m, y_m, h_m, tick);
        let mut velocity = turb.u;
        velocity[0] += self.mean.speed(h_m);
        let mut grad = turb.grad;
        grad[0][2] += self.mean.dspeed_dh(h_m);
        Ok(AirState {
            velocity_mps: velocity,
            grad,
            rho_kg_m3: self.air.rho_kg_m3,
            mu_kg_m_s: self.air.mu_kg_m_s,
            temperature_k: self.air.temperature_k,
            pressure_pa: self.air.pressure_pa,
            provenance: self.air.provenance,
        })
    }
}
