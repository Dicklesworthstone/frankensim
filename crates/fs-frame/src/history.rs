//! Stage 3: nonlinear TIME HISTORY of a single-story frame with two
//! fiber-hinge columns — the smoke-tier concentrated-plasticity
//! idealization: story drift x maps to base-hinge curvature
//! κ = x/(h·l_p), the TRUE fiber section (Mander concrete core +
//! Menegotto–Pinto steel through fs-solid/fs-material, with all the
//! sign conventions tfz.14 pinned) returns the hinge moment, and the
//! story shear is V = 2M/h. Newmark average acceleration with Newton
//! on the section tangent; hysteresis and stiffness degradation come
//! from the fibers, not from a phenomenological spring. The
//! distributed-plasticity frame (fs-solid ForceBasedElement columns)
//! is the recorded successor.

use fs_solid::fiber::{Section, rc_section};

/// Story model parameters.
#[derive(Debug, Clone, Copy)]
pub struct StoryParams {
    /// Story height (m).
    pub h: f64,
    /// Plastic-hinge length (m).
    pub lp: f64,
    /// Story mass (t · consistent units).
    pub mass: f64,
    /// Damping ratio (Rayleigh mass-proportional at the initial
    /// period).
    pub zeta: f64,
    /// Section scale: fiber areas multiply by this (the CVaR design
    /// variable).
    pub scale: f64,
}

impl Default for StoryParams {
    fn default() -> Self {
        StoryParams {
            h: 3.0,
            lp: 0.45,
            // 280 t story mass over the two-column pair: T ≈ 0.5 s at
            // the probed k₀ ≈ 4.4e7 N/m; yield drift ratio ≈ 0.0036
            // (V_y ≈ 4.8e5 N). SI throughout — fs-material's units.
            mass: 2.8e5,
            zeta: 0.02,
            scale: 1.0,
        }
    }
}

/// The story frame state (two identical fiber-hinge columns).
pub struct StoryFrame {
    /// Parameters.
    pub params: StoryParams,
    hinge: Section,
    /// Committed drift.
    pub x: f64,
    /// Committed velocity.
    pub v: f64,
    /// Committed acceleration (relative).
    pub a: f64,
}

/// Scale a section's fiber areas.
fn scaled_section(scale: f64) -> Section {
    let mut s = rc_section(0.5, 0.35, 12, 0.002);
    for f in &mut s.fibers {
        f.area *= scale;
    }
    s
}

impl StoryFrame {
    /// A story at rest.
    #[must_use]
    pub fn new(params: StoryParams) -> StoryFrame {
        StoryFrame {
            params,
            hinge: scaled_section(params.scale),
            x: 0.0,
            v: 0.0,
            a: 0.0,
        }
    }

    /// Story restoring shear and tangent stiffness at drift `x`
    /// (trial — no commit).
    #[must_use]
    pub fn restoring(&self, x: f64) -> (f64, f64) {
        let kappa = x / (self.params.h * self.params.lp);
        let st = self.hinge.respond(0.0, kappa);
        let v = 2.0 * st.m / self.params.h;
        let dv_dx = 2.0 * st.tangent[1][1] / (self.params.h * self.params.h * self.params.lp);
        (v, dv_dx)
    }

    /// Initial (elastic) story stiffness.
    #[must_use]
    pub fn initial_stiffness(&self) -> f64 {
        self.restoring(1e-9).1
    }

    /// Run the record under ground acceleration `ag` sampled at `dt`;
    /// returns the drift history. Newmark average acceleration
    /// (γ = ½, β = ¼) with Newton on the fiber tangent; the section
    /// commits once per step.
    pub fn run(&mut self, ag: &[f64], dt: f64) -> Vec<f64> {
        let m = self.params.mass;
        let k0 = self.initial_stiffness();
        let c = 2.0 * self.params.zeta * (k0 * m).sqrt();
        let (beta, gamma) = (0.25f64, 0.5f64);
        let mut drifts = Vec::with_capacity(ag.len());
        for &agi in ag {
            let p_ext = -m * agi;
            // Newmark predictors.
            let x0 = self.x;
            let v0 = self.v;
            let a0 = self.a;
            let mut x = x0;
            for _ in 0..30 {
                let a_new = (x - x0 - dt * v0) / (beta * dt * dt) - (0.5 - beta) / beta * a0;
                let v_new = v0 + dt * ((1.0 - gamma) * a0 + gamma * a_new);
                let (fs, kt) = self.restoring(x);
                let r = m * a_new + c * v_new + fs - p_ext;
                let kdyn = m / (beta * dt * dt) + c * gamma / (beta * dt) + kt;
                let dx = -r / kdyn;
                x += dx;
                if dx.abs() < 1e-12 {
                    break;
                }
            }
            let a_new = (x - x0 - dt * v0) / (beta * dt * dt) - (0.5 - beta) / beta * a0;
            let v_new = v0 + dt * ((1.0 - gamma) * a0 + gamma * a_new);
            let kappa = x / (self.params.h * self.params.lp);
            self.hinge.commit(0.0, kappa);
            self.x = x;
            self.v = v_new;
            self.a = a_new;
            drifts.push(x);
        }
        drifts
    }
}

/// Peak drift RATIO (|x|max / h) of a drift history.
#[must_use]
pub fn peak_drift(drifts: &[f64], h: f64) -> f64 {
    drifts.iter().fold(0.0f64, |m, &x| m.max(x.abs())) / h
}
