//! Normal-flow, nonrigid-wake first harmonics of actuator-disk inflow.
//!
//! For `v_i(x,y) = u_0 + h_x x/R + h_y y/R`, the conjugate blade loads are
//! `[T, -M_y/R, M_x/R]` in a right-handed disk frame. The dimensional masses
//! are `m_0 = (8/3) rho R^3`, `m_1 = (16/45) rho R^3`. Normal mean transport
//! is the original Carpenter--Fridovich law; harmonic transport is
//! `rho A (V+2u_0) h / 2`. These are the zero-skew, nonrigid-wake limits of
//! Chen, NASA TM-102219, equations (13b), (14), (17), (18) and (21):
//! <https://ntrs.nasa.gov/api/citations/19900006622/downloads/19900006622.pdf>.
//!
//! Both harmonics use a nonspinning material disk basis. Its rotation about
//! the disk normal adds the coordinate term `m_1 omega_n J h`, with
//! `J[h_x,h_y]=[-h_y,h_x]`. This term does no work; rotor shaft speed is NOT
//! the basis rate. Tangent transport follows the current disk plane, not a
//! resolved three-dimensional wake. `skew` separately selects the full skew
//! transport. Wake distortion, ground effect and higher spatial harmonics
//! remain outside both choices.

use super::{Refusal, UniformAxialInflow, invalid};

/// Full first-harmonic skew transport in fixed disk coordinates.
pub mod skew;

/// Dimensional mean and tip-amplitude spatial modes, all in m/s.
/// The ordered modes are mean, cosine along the disk radial reference, and
/// sine along `axis cross reference`. They do not rotate with blade phase.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicInduction {
    pub mean_mps: f64,
    pub first_mps: [f64; 2],
}

impl HarmonicInduction {
    #[must_use]
    pub fn modes(self) -> [f64; 3] {
        [self.mean_mps, self.first_mps[0], self.first_mps[1]]
    }

    /// Point evaluation; `x_over_r` and `y_over_r` are coordinates divided by
    /// the outer disk radius, not by the local radial station.
    #[must_use]
    pub fn at(self, x_over_r: f64, y_over_r: f64) -> f64 {
        self.mean_mps + self.first_mps[0] * x_over_r + self.first_mps[1] * y_over_r
    }
}

/// Constant-density, three-mode normal-flow model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FirstHarmonicInflow {
    uniform: UniformAxialInflow,
    mass_kg: [f64; 3],
}

/// Continuous generalized force balance, each force in N.
/// Harmonic forces are transverse moment divided by disk radius.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicInflowRate {
    pub acceleration_mps2: [f64; 3],
    pub transport_force_n: [f64; 3],
    /// Moving-coordinate term; not a physical source or aerodynamic load.
    pub frame_force_n: [f64; 3],
}

/// Implicit-midpoint modal step under held blade loads, incident velocity
/// and disk-basis axial rate. These are reduced modal balances,
/// not the world momentum or total kinetic energy of a resolved wake.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicInflowStep {
    pub induction: HarmonicInduction,
    pub midpoint: HarmonicInduction,
    pub driving_impulse_ns: [f64; 3],
    pub transport_impulse_ns: [f64; 3],
    pub frame_impulse_ns: [f64; 3],
    /// `M delta_u - driving + transport + frame`.
    pub momentum_residual_ns: [f64; 3],
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Generalized blade force dotted with midpoint modal velocity.
    pub induced_work_j: f64,
    /// Actual transport work; skew mean/harmonic coupling has a signed part.
    pub transport_work_j: f64,
    /// Zero in real arithmetic, retained to expose floating-point transport.
    pub frame_work_j: f64,
    pub energy_residual_j: f64,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]
}

impl FirstHarmonicInflow {
    /// Prepare the potential-flow modal masses. Neither mass is added to
    /// the solid rotor or carrier inertia.
    pub fn new(radius_m: f64, density_kg_m3: f64) -> Result<Self, Refusal> {
        let uniform = UniformAxialInflow::new(radius_m, density_kg_m3)?;
        let first_mass = (16.0 / 45.0) * density_kg_m3 * radius_m * radius_m * radius_m;
        if !first_mass.is_finite() || first_mass <= 0.0 {
            return Err(invalid("unrepresentable harmonic apparent mass"));
        }
        Ok(Self {
            mass_kg: [uniform.apparent_mass_kg(), first_mass, first_mass],
            uniform,
        })
    }

    #[must_use]
    pub fn density_kg_m3(self) -> f64 {
        self.uniform.density_kg_m3()
    }

    #[must_use]
    pub fn modal_mass_kg(self) -> [f64; 3] {
        self.mass_kg
    }

    fn admit(self, state: HarmonicInduction, axial_mps: f64) -> Result<(), Refusal> {
        self.uniform.admit_flow(axial_mps, state.mean_mps)?;
        if !state.first_mps.into_iter().all(f64::is_finite)
            || state.first_mps[0].hypot(state.first_mps[1]) > axial_mps + state.mean_mps
        {
            return Err(invalid(
                "harmonic induction reverses normal throughflow somewhere on the disk",
            ));
        }
        Ok(())
    }

    fn damping(self, axial_mps: f64, mean_mps: f64) -> f64 {
        0.5 * self.uniform.density_kg_m3 * self.uniform.area_m2 * (axial_mps + 2.0 * mean_mps)
    }

    /// Evaluate actual modal acceleration without a steady wake solve.
    /// `frame_rate_rad_s` is carrier rotation about the disk normal. It
    /// excludes relative rotor rotation and does not create physical work.
    pub fn rate(
        self,
        state: HarmonicInduction,
        axial_mps: f64,
        force_n: [f64; 3],
        frame_rate_rad_s: f64,
    ) -> Result<HarmonicInflowRate, Refusal> {
        self.admit(state, axial_mps)?;
        let mean = self.uniform.rate(axial_mps, state.mean_mps, force_n[0])?;
        let damping = self.damping(axial_mps, state.mean_mps);
        let transport_force_n = [
            mean.momentum_flux_n,
            damping * state.first_mps[0],
            damping * state.first_mps[1],
        ];
        let frame_force_n = [
            0.0,
            -self.mass_kg[1] * frame_rate_rad_s * state.first_mps[1],
            self.mass_kg[1] * frame_rate_rad_s * state.first_mps[0],
        ];
        let acceleration_mps2 = core::array::from_fn(|j| {
            (force_n[j] - transport_force_n[j] - frame_force_n[j]) / self.mass_kg[j]
        });
        if !force_n
            .iter()
            .chain(&transport_force_n)
            .chain(&frame_force_n)
            .chain(&acceleration_mps2)
            .chain(&[frame_rate_rad_s, damping])
            .all(|x| x.is_finite())
        {
            return Err(invalid("harmonic force or rate overflow"));
        }
        Ok(HarmonicInflowRate {
            acceleration_mps2,
            transport_force_n,
            frame_force_n,
        })
    }

    /// Solve the nonlinear mean quadratic and both coupled harmonic midpoint
    /// equations analytically. Constant work/storage; no iteration or allocation.
    /// Actual initial, midpoint and final whole-disk throughflow must be admitted.
    /// Midpoint is second order and dissipative for held normal-flow data, but
    /// is not L-stable and does not resolve an arbitrarily coarse transient.
    pub fn step_held(
        self,
        before: HarmonicInduction,
        axial_mps: f64,
        force_n: [f64; 3],
        frame_rate_rad_s: f64,
        dt_s: f64,
    ) -> Result<HarmonicInflowStep, Refusal> {
        self.rate(before, axial_mps, force_n, frame_rate_rad_s)?;
        let mean = self
            .uniform
            .step_held(before.mean_mps, axial_mps, force_n[0], dt_s)?;
        let a = 0.5 * dt_s * self.damping(axial_mps, mean.midpoint_induced_mps) / self.mass_kg[1];
        let b = 0.5 * dt_s * frame_rate_rad_s;
        let c = 1.0 + a;
        let denominator = c * c + b * b;
        let rhs = [
            (1.0 - a) * before.first_mps[0]
                + b * before.first_mps[1]
                + dt_s * force_n[1] / self.mass_kg[1],
            (1.0 - a) * before.first_mps[1] - b * before.first_mps[0]
                + dt_s * force_n[2] / self.mass_kg[1],
        ];
        if ![a, b, c, denominator, rhs[0], rhs[1]]
            .into_iter()
            .all(f64::is_finite)
            || denominator <= 0.0
        {
            return Err(invalid("unrepresentable harmonic midpoint equation"));
        }
        let induction = HarmonicInduction {
            mean_mps: mean.induced_mps,
            first_mps: [
                (c * rhs[0] + b * rhs[1]) / denominator,
                (c * rhs[1] - b * rhs[0]) / denominator,
            ],
        };
        let midpoint = HarmonicInduction {
            mean_mps: mean.midpoint_induced_mps,
            first_mps: core::array::from_fn(|j| {
                0.5 * (before.first_mps[j] + induction.first_mps[j])
            }),
        };
        self.admit(induction, axial_mps)?;
        let rate = self.rate(midpoint, axial_mps, force_n, frame_rate_rad_s)?;
        self.account_step(before, induction, midpoint, force_n, rate, dt_s)
    }

    // Both physical transport choices use this same original storage/impulse owner.
    fn account_step(
        self,
        before: HarmonicInduction,
        induction: HarmonicInduction,
        midpoint: HarmonicInduction,
        force_n: [f64; 3],
        rate: HarmonicInflowRate,
        dt_s: f64,
    ) -> Result<HarmonicInflowStep, Refusal> {
        let driving_impulse_ns = force_n.map(|f| dt_s * f);
        let transport_impulse_ns = rate.transport_force_n.map(|f| dt_s * f);
        let frame_impulse_ns = rate.frame_force_n.map(|f| dt_s * f);
        let delta = core::array::from_fn::<_, 3, _>(|j| induction.modes()[j] - before.modes()[j]);
        let momentum_residual_ns = core::array::from_fn(|j| {
            self.mass_kg[j] * delta[j] - driving_impulse_ns[j]
                + transport_impulse_ns[j]
                + frame_impulse_ns[j]
        });
        let kinetic =
            |state: HarmonicInduction| 0.5 * dot(self.mass_kg, state.modes().map(|x| x * x));
        let kinetic_before_j = kinetic(before);
        let kinetic_after_j = kinetic(induction);
        let induced_work_j = dot(driving_impulse_ns, midpoint.modes());
        let transport_work_j = dot(transport_impulse_ns, midpoint.modes());
        let frame_work_j = dot(frame_impulse_ns, midpoint.modes());
        let energy_residual_j =
            kinetic_after_j - kinetic_before_j - induced_work_j + transport_work_j + frame_work_j;
        if !driving_impulse_ns
            .iter()
            .chain(&transport_impulse_ns)
            .chain(&frame_impulse_ns)
            .chain(&momentum_residual_ns)
            .chain(&[
                kinetic_before_j,
                kinetic_after_j,
                induced_work_j,
                transport_work_j,
                frame_work_j,
                energy_residual_j,
            ])
            .all(|x| x.is_finite())
        {
            return Err(invalid("harmonic midpoint accounting overflow"));
        }
        Ok(HarmonicInflowStep {
            induction,
            midpoint,
            driving_impulse_ns,
            transport_impulse_ns,
            frame_impulse_ns,
            momentum_residual_ns,
            kinetic_before_j,
            kinetic_after_j,
            induced_work_j,
            transport_work_j,
            frame_work_j,
            energy_residual_j,
        })
    }
}
