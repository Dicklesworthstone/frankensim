//! Uniform axial dynamic inflow with an explicit apparent-air-mass state.
//!
//! Carpenter–Fridovich's normal-flow model is
//! `m_a dv_i/dt = T - 2 rho A v_i (V + v_i)`, with
//! `m_a = (8/3) rho R^3`. See Chen and Hindson, NASA TM-88327,
//! equation (1), <https://ntrs.nasa.gov/citations/19870015049>.
//! The coefficient is the potential-flow apparent mass of a normal disk,
//! not a fitted time constant or additional rigid-aircraft mass.
//! Positive `v_i` increases the oncoming relative flow seen by the blades;
//! the corresponding physical air impulse opposes aircraft thrust.
//!
//! This is a declared UNIFORM-disk model. Blade forces reuse the native
//! section law and radial quadrature at the current induced velocity;
//! they are never obtained by solving a steady BEMT problem and delaying
//! its answer. Uniform disk momentum replaces the steady annular Prandtl
//! closure, so their steady radial distributions are different. No blade
//! unsteadiness, azimuthal harmonics, skew wake, swirl, reverse throughflow,
//! moving-wake geometry or historical validation is supplied.

use super::{MAX_BRAKING_INDUCTION, Refusal, Rotor, StationProblem, refuse};

/// Mean and first spatial harmonics with physical modal air inertia.
pub mod first_harmonic;

/// Constant-density uniform disk, with SI geometry and physical air inertia.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniformAxialInflow {
    radius_m: f64,
    density_kg_m3: f64,
    area_m2: f64,
    apparent_mass_kg: f64,
}

/// Continuous axial momentum balance at one stored induced velocity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InflowRate {
    /// `dv_i/dt` [m/s²].
    pub acceleration_mps2: f64,
    /// Advected axial momentum per unit time `2 rho A v_i (V+v_i)` [N].
    pub momentum_flux_n: f64,
}

/// Actual section loads and wake acceleration at the same physical state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorEvaluation {
    /// Thrust along the positive disk axis [N], without sign clipping.
    pub thrust_n: f64,
    /// Required torque in the positive blade-spin convention [N m].
    pub torque_nm: f64,
    /// Wake momentum balance; a nonzero rate is a physical transient.
    pub inflow: InflowRate,
}

/// Implicit-midpoint wake step under held thrust and incident axial speed.
/// Work here is thrust times INDUCED velocity; it is not total shaft work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InflowStep {
    /// Accepted disk-induced velocity [m/s].
    pub induced_mps: f64,
    /// Arithmetic midpoint used by the implicit wake flux [m/s].
    pub midpoint_induced_mps: f64,
    /// Blade impulse into the induced-flow coordinate [N s]; the physical
    /// air reaction is opposite the positive aircraft-thrust direction.
    pub thrust_impulse_ns: f64,
    /// Momentum transported out of this reduced wake state [N s].
    pub transport_impulse_ns: f64,
    /// `m_a Δv_i - thrust_impulse + transport_impulse` [N s].
    pub momentum_residual_ns: f64,
    /// Apparent-fluid kinetic energy before this step [J].
    pub kinetic_before_j: f64,
    /// Apparent-fluid kinetic energy after this step [J].
    pub kinetic_after_j: f64,
    /// Held blade thrust times midpoint induced velocity times duration [J].
    pub induced_work_j: f64,
    /// Work carried by the modeled momentum flux [J].
    pub transport_work_j: f64,
    /// Change in stored energy minus induced work plus transport work [J].
    pub energy_residual_j: f64,
}

fn invalid(message: &str) -> Refusal {
    refuse(
        "dynamic-inflow-invalid",
        message.into(),
        "supply finite constant-density normal axial flow and a matching physical rotor",
    )
}

impl UniformAxialInflow {
    /// Prepare the physical disk area and its apparent air mass.
    ///
    /// # Errors
    /// `dynamic-inflow-invalid` for nonpositive/nonfinite or overflowing data.
    pub fn new(radius_m: f64, density_kg_m3: f64) -> Result<Self, Refusal> {
        let area_m2 = core::f64::consts::PI * radius_m * radius_m;
        let apparent_mass_kg = (8.0 / 3.0) * density_kg_m3 * radius_m * radius_m * radius_m;
        if ![radius_m, density_kg_m3, area_m2, apparent_mass_kg]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
        {
            return Err(invalid("invalid radius, density, area or apparent mass"));
        }
        Ok(Self {
            radius_m,
            density_kg_m3,
            area_m2,
            apparent_mass_kg,
        })
    }

    /// Fixed density [kg/m³]. It must not change while this wake is advanced.
    #[must_use]
    pub fn density_kg_m3(self) -> f64 {
        self.density_kg_m3
    }

    /// Apparent normal air mass [kg], excluded from carrier mass/inertia.
    #[must_use]
    pub fn apparent_mass_kg(self) -> f64 {
        self.apparent_mass_kg
    }

    fn admit_flow(self, axial_mps: f64, induced_mps: f64) -> Result<(), Refusal> {
        if !axial_mps.is_finite() || !induced_mps.is_finite() || axial_mps < 0.0 {
            return Err(invalid(
                "nonfinite flow or reversed incident axial velocity",
            ));
        }
        if induced_mps < -MAX_BRAKING_INDUCTION * axial_mps {
            return Err(refuse(
                "dynamic-inflow-wake-domain",
                "wake is outside the normal-flow branch (braking induction exceeds 0.4)".into(),
                "reduce the timestep if the transient overshot; otherwise supply an applicable wake model",
            ));
        }
        Ok(())
    }

    /// Evaluate the continuous wake momentum equation for a declared thrust.
    /// Signed thrust and moderate braking are retained without clamping.
    ///
    /// # Errors
    /// `dynamic-inflow-invalid` or `dynamic-inflow-wake-domain`.
    pub fn rate(
        self,
        axial_mps: f64,
        induced_mps: f64,
        thrust_n: f64,
    ) -> Result<InflowRate, Refusal> {
        self.admit_flow(axial_mps, induced_mps)?;
        let momentum_flux_n =
            2.0 * self.density_kg_m3 * self.area_m2 * induced_mps * (axial_mps + induced_mps);
        let acceleration_mps2 = (thrust_n - momentum_flux_n) / self.apparent_mass_kg;
        if ![thrust_n, momentum_flux_n, acceleration_mps2]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("unrepresentable wake momentum equation"));
        }
        Ok(InflowRate {
            acceleration_mps2,
            momentum_flux_n,
        })
    }

    /// Integrate the SAME blade section force law used by `bemt_solve` at a
    /// stored, spatially uniform induction. No annular equilibrium is imposed.
    ///
    /// # Errors
    /// Native rotor/section admission, `dynamic-inflow-invalid`, or
    /// `dynamic-inflow-wake-domain`. Positive physical blade speed is required.
    pub fn evaluate(
        self,
        rotor: &Rotor,
        axial_mps: f64,
        induced_mps: f64,
        omega_rad_s: f64,
    ) -> Result<RotorEvaluation, Refusal> {
        rotor.admit()?;
        self.admit_flow(axial_mps, induced_mps)?;
        if rotor.radius_m != self.radius_m || !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid("rotor radius mismatch or nonpositive blade speed"));
        }
        let root_radius = rotor.stations[0].r_over_r * 0.8 * rotor.radius_m;
        let mut thrust_n = 0.0;
        let mut torque_nm = 0.0;
        for (index, station) in rotor.stations.iter().enumerate() {
            let radius = station.r_over_r * rotor.radius_m;
            let loads = StationProblem {
                rotor,
                station,
                index,
                rho: self.density_kg_m3,
                axial: axial_mps,
                tangential: omega_rad_s * radius,
                radius,
                root_radius,
            }
            .blade_loads(induced_mps)?;
            let width = rotor.station_width_m(index);
            thrust_n += loads.thrust * width;
            torque_nm += loads.torque * width;
        }
        if !torque_nm.is_finite() {
            return Err(invalid("unrepresentable integrated shaft torque"));
        }
        Ok(RotorEvaluation {
            thrust_n,
            torque_nm,
            inflow: self.rate(axial_mps, induced_mps, thrust_n)?,
        })
    }

    /// Advance the nonlinear wake flux with implicit midpoint. The scalar
    /// quadratic is solved analytically, without an iterative wake solver.
    /// Its rationalized normal root retains arbitrarily weak loading.
    /// Momentum and the quadratic storage/work balance close to roundoff.
    /// This is second order for held inputs, not a global coupled passivity
    /// claim. A step that leaves normal throughflow refuses instead of clipping.
    ///
    /// # Errors
    /// `dynamic-inflow-invalid`, `dynamic-inflow-wake-domain`, or
    /// `dynamic-inflow-step-domain` if the midpoint quadratic has no real root.
    pub fn step_held(
        self,
        induced_mps: f64,
        axial_mps: f64,
        thrust_n: f64,
        dt_s: f64,
    ) -> Result<InflowStep, Refusal> {
        self.rate(axial_mps, induced_mps, thrust_n)?;
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid("wake step duration must be positive finite"));
        }
        let half_step_per_mass = 0.5 * dt_s / self.apparent_mass_kg;
        let a = half_step_per_mass * 2.0 * self.density_kg_m3 * self.area_m2;
        let b = 1.0 + a * axial_mps;
        let c = induced_mps + half_step_per_mass * thrust_n;
        let discriminant = b * b + 4.0 * a * c;
        if !discriminant.is_finite() || discriminant < 0.0 {
            return Err(refuse(
                "dynamic-inflow-step-domain",
                "normal midpoint wake root is non-real or unrepresentable".into(),
                "reduce the timestep or use an applicable non-normal-flow wake model",
            ));
        }
        let midpoint = 2.0 * c / (b + discriminant.sqrt());
        let after = midpoint + (midpoint - induced_mps);
        self.admit_flow(axial_mps, after)?;
        let flux = self.rate(axial_mps, midpoint, thrust_n)?.momentum_flux_n;
        let thrust_impulse_ns = dt_s * thrust_n;
        let transport_impulse_ns = dt_s * flux;
        let momentum_residual_ns = self.apparent_mass_kg * (after - induced_mps)
            - thrust_impulse_ns
            + transport_impulse_ns;
        let kinetic_before_j = 0.5 * self.apparent_mass_kg * induced_mps * induced_mps;
        let kinetic_after_j = 0.5 * self.apparent_mass_kg * after * after;
        let induced_work_j = thrust_impulse_ns * midpoint;
        let transport_work_j = transport_impulse_ns * midpoint;
        let energy_residual_j =
            kinetic_after_j - kinetic_before_j - induced_work_j + transport_work_j;
        if ![
            midpoint,
            after,
            thrust_impulse_ns,
            transport_impulse_ns,
            momentum_residual_ns,
            kinetic_before_j,
            kinetic_after_j,
            induced_work_j,
            transport_work_j,
            energy_residual_j,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(invalid("unrepresentable wake impulse or work"));
        }
        Ok(InflowStep {
            induced_mps: after,
            midpoint_induced_mps: midpoint,
            thrust_impulse_ns,
            transport_impulse_ns,
            momentum_residual_ns,
            kinetic_before_j,
            kinetic_after_j,
            induced_work_j,
            transport_work_j,
            energy_residual_j,
        })
    }
}
