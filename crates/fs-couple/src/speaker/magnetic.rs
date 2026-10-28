//! Excursion-dependent voice-coil magnetics from one reversible energy.
//!
//! Total linkage lambda is the electrical state, with
//! H_m = (lambda - psi_pm(x))^2 / (2 L(x)). Thus current is dH/dlambda and
//! force is -dH/dx = Bl(x) i + L'(x) i^2/2. Since psi_pm' = Bl,
//! lambda_dot = voltage - R i includes BOTH permanent-field back EMF and
//! motional inductance voltage. A separately patched force or voltage would
//! break this reciprocity. Coefficients are caller data, not calibrated claims.
//! The reversible Bl(x)/Le(x) terms follow the electrodynamic model in
//! Klippel, "Predicting the Nonlinear Loudspeaker Behavior", Training 5:
//! <https://www.klippel.de/training/attachments/training5/Training_5_Predicting_the_Nonlinear_Loudspeaker_Behavior_en.pdf>.
//! Its other current-dependent, thermal and lossy motor terms are not supplied
//! by this displacement-only profile.

/// Source-bound voice-coil heating on the shared finite enthalpy timeline.
pub mod thermal;

use super::{DriverAssembly, PanelMode, SpeakerError, TsCard, assemble_driver, declared_air};
use fs_material::gas::GasState;
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage, StepRecord, Storage, step_descriptor};

fn invalid(what: &'static str) -> SpeakerError {
    SpeakerError::Invalid { what }
}

/// Quadratic Bl(x) and L(x) profiles on an explicitly bounded displacement range.
///
/// With xi = x / displacement_scale_m, each coefficient array means
/// c0 + c1 xi + c2 xi^2. All Bl coefficients have units [T m]; all L
/// coefficients have units [H]. These local measured/field-fitted profiles
/// retain asymmetry and curvature without inventing extrapolation or data.
/// Permanent linkage is their analytic integral with psi_pm(0) = 0.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VoiceCoilMagnetics {
    displacement_limits_m: [f64; 2],
    displacement_scale_m: f64,
    force_factor_tm: [f64; 3],
    inductance_h: [f64; 3],
}

/// One admitted local magnetic constitutive evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MagneticCoilSample {
    /// Signed permanent-magnet linkage relative to x=0 [Wb turns].
    pub permanent_linkage_wb: f64,
    /// Reciprocal permanent-field force/back-EMF factor [T m].
    pub force_factor_tm: f64,
    /// Positive differential self-inductance [H].
    pub inductance_h: f64,
    /// Spatial derivative of inductance [H/m], which supplies reluctance force.
    pub inductance_gradient_h_per_m: f64,
}

impl VoiceCoilMagnetics {
    /// Admit finite coefficients, a positive length scale, a range containing
    /// zero, and strictly positive L over the WHOLE interval (including its
    /// interior minimum). No clipping of a negative-inductance fit is allowed.
    pub fn new(
        displacement_limits_m: [f64; 2],
        displacement_scale_m: f64,
        force_factor_tm: [f64; 3],
        inductance_h: [f64; 3],
    ) -> Result<Self, SpeakerError> {
        let [lo, hi] = displacement_limits_m;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo >= 0.0
            || hi <= 0.0
            || !displacement_scale_m.is_finite()
            || displacement_scale_m <= 0.0
            || force_factor_tm
                .iter()
                .chain(&inductance_h)
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "magnetic profile needs finite data and a displacement range containing zero",
            ));
        }
        let profile = Self {
            displacement_limits_m,
            displacement_scale_m,
            force_factor_tm,
            inductance_h,
        };
        let extent = lo.abs().max(hi.abs()) / displacement_scale_m;
        // Finite absolute term bounds rule out hidden intermediate overflows
        // inside a polynomial whose endpoints happen to cancel.
        for c in [force_factor_tm, inductance_h] {
            if !(c[0].abs() + extent * (c[1].abs() + extent * c[2].abs())).is_finite()
                || !((c[1].abs() + 2.0 * extent * c[2].abs()) / displacement_scale_m).is_finite()
            {
                return Err(invalid(
                    "magnetic profile is not representable over its displacement range",
                ));
            }
        }
        for x in [lo, 0.0, hi] {
            profile.sample(x)?;
        }
        if inductance_h[2] > 0.0 {
            let vertex = -0.5 * inductance_h[1] / inductance_h[2] * displacement_scale_m;
            if vertex > lo && vertex < hi {
                profile.sample(vertex)?;
            }
        }
        Ok(profile)
    }

    /// Closed validity interval [m]. No value outside it is extrapolated.
    pub const fn displacement_limits_m(&self) -> [f64; 2] {
        self.displacement_limits_m
    }

    /// Evaluate the admitted constitutive law at a physical displacement [m].
    pub fn sample(&self, x: f64) -> Result<MagneticCoilSample, SpeakerError> {
        if !x.is_finite() || x < self.displacement_limits_m[0] || x > self.displacement_limits_m[1]
        {
            return Err(invalid("coil displacement is outside the magnetic profile"));
        }
        let xi = x / self.displacement_scale_m;
        let b = self.force_factor_tm;
        let l = self.inductance_h;
        let sample = MagneticCoilSample {
            permanent_linkage_wb: x * (b[0] + xi * (b[1] / 2.0 + xi * b[2] / 3.0)),
            force_factor_tm: b[0] + xi * (b[1] + xi * b[2]),
            inductance_h: l[0] + xi * (l[1] + xi * l[2]),
            inductance_gradient_h_per_m: (l[1] + 2.0 * xi * l[2]) / self.displacement_scale_m,
        };
        if !sample.permanent_linkage_wb.is_finite()
            || !sample.force_factor_tm.is_finite()
            || !sample.inductance_h.is_finite()
            || sample.inductance_h <= 0.0
            || !sample.inductance_gradient_h_per_m.is_finite()
        {
            return Err(invalid(
                "magnetic profile requires finite linkage and strictly positive inductance",
            ));
        }
        Ok(sample)
    }
}

#[derive(Clone)]
struct CoilStorage {
    mechanical: QuadraticStorage,
    profile: VoiceCoilMagnetics,
    exact_q: Option<Vec<f64>>,
}

impl Storage for CoilStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        let Ok(local) = self.profile.sample(x[1]) else {
            return f64::NAN;
        };
        let flux = x[0] - local.permanent_linkage_wb;
        self.mechanical.hamiltonian(x) + 0.5 * flux * (flux / local.inductance_h)
    }

    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        let Ok(local) = self.profile.sample(x[1]) else {
            out.fill(f64::NAN);
            return;
        };
        self.mechanical.gradient(x, out);
        let current = (x[0] - local.permanent_linkage_wb) / local.inductance_h;
        out[0] += current;
        out[1] -= current * local.force_factor_tm
            + 0.5 * current * current * local.inductance_gradient_h_per_m;
    }

    fn exact_quadratic_matrix(&self) -> Option<&[f64]> {
        self.exact_q.as_deref()
    }

    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        self.exact_q
            .as_ref()
            .map(|q| fs_phs::linear_exact::QuadraticPart {
                matrix: q.clone(),
                exact: true,
            })
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let (Ok(la), Ok(lb), Ok(lm)) = (
            self.profile.sample(a[1]),
            self.profile.sample(b[1]),
            self.profile.sample(mid[1]),
        ) else {
            return Some(f64::NAN);
        };
        let ia = (a[0] - la.permanent_linkage_wb) / la.inductance_h;
        let ib = (b[0] - lb.permanent_linkage_wb) / lb.inductance_h;
        let im = (mid[0] - lm.permanent_linkage_wb) / lm.inductance_h;
        let ua = a[1] / self.profile.displacement_scale_m;
        let ub = b[1] / self.profile.displacement_scale_m;
        let c = self.profile.force_factor_tm;
        let psi_secant = c[0] + 0.5 * c[1] * (ua + ub) + c[2] / 3.0 * (ua * ua + ua * ub + ub * ub);
        let l = self.profile.inductance_h;
        let l_secant = (l[1] + l[2] * (ua + ub)) / self.profile.displacement_scale_m;
        // Delta Hm = (ia+ib)/2 * Delta(lambda-psi) - ia*ib/2 * Delta L.
        // Subtract the midpoint gradient algebraically; no large endpoint
        // energies or nearly equal permanent linkage values are subtracted.
        let avg_i = f64::midpoint(ia, ib);
        Some(
            (avg_i - im) * (b[0] - a[0])
                + (im * lm.force_factor_tm + 0.5 * im * im * lm.inductance_gradient_h_per_m
                    - avg_i * psi_secant
                    - 0.5 * ia * ib * l_secant)
                    * (b[1] - a[1]),
        )
    }
}

/// Opt-in nonlinear magnetic speaker, retaining the original TS cone, sealed
/// cavity and optional panel mechanics. State 0 is TOTAL linkage, not L*i.
///
/// The existing descriptor discrete-gradient solver advances electrical and
/// mechanical states jointly. It conserves the original discrete energy/work
/// balance; time accuracy remains second order, unlike the linear TsDriver's
/// exact held-input propagator. Resolve the nonlinear response by timestep
/// refinement. This small dense reference owner is not an audio callback.
pub struct MagneticTsDriver {
    phs: DescriptorPortHamiltonian,
    storage: CoilStorage,
    profile: VoiceCoilMagnetics,
    x: Vec<f64>,
    card: TsCard,
    rho: f64,
    moving_mass_kg: f64,
    has_panel: bool,
    time_s: f64,
}

impl MagneticTsDriver {
    /// Compose a nonlinear magnetic driver in the same declared room air as
    /// [`super::TsDriver::new`]. The profile's origin must match card Bl/Le.
    pub fn new(
        card: TsCard,
        box_volume_m3: Option<f64>,
        panel: Option<PanelMode>,
        profile: VoiceCoilMagnetics,
    ) -> Result<Self, SpeakerError> {
        Self::new_with_gas(card, box_volume_m3, panel, declared_air(), profile)
    }

    /// Compose with explicit gas properties. The magnetic validity interval
    /// must fit inside the card's linear suspension excursion. No suspension,
    /// core saturation, hysteresis, skin loss or coil-heating law is inferred.
    pub fn new_with_gas(
        card: TsCard,
        box_volume_m3: Option<f64>,
        panel: Option<PanelMode>,
        gas: GasState,
        profile: VoiceCoilMagnetics,
    ) -> Result<Self, SpeakerError> {
        let DriverAssembly {
            n,
            mut q,
            mut j,
            r,
            g,
            moving_mass_kg,
        } = assemble_driver(card, box_volume_m3, panel, &gas)?;
        let origin = profile.sample(0.0)?;
        if (origin.force_factor_tm - card.bl).abs() > 1e-12 * card.bl
            || (origin.inductance_h - card.le_h).abs() > 1e-12 * card.le_h
            || profile.displacement_limits_m[0] < -card.x_max_m
            || profile.displacement_limits_m[1] > card.x_max_m
        {
            return Err(invalid(
                "magnetic profile must match card Bl/Le at zero and stay within suspension travel",
            ));
        }
        // Total-linkage coordinates put ALL transduction into H. Retaining
        // the linear Bl gyrator too would double both force and back EMF.
        j[2] = 0.0;
        j[2 * n] = 0.0;
        q[0] = 0.0;
        let exact_q = if profile.force_factor_tm[1..] == [0.0, 0.0]
            && profile.inductance_h[1..] == [0.0, 0.0]
        {
            let mut full = q.clone();
            full[0] = 1.0 / card.le_h;
            full[1] = -card.bl / card.le_h;
            full[n] = full[1];
            full[n + 1] += card.bl * card.bl / card.le_h;
            Some(full)
        } else {
            None
        };
        let storage = CoilStorage {
            mechanical: QuadraticStorage::new(q, n)?,
            profile,
            exact_q,
        };
        let phs = DescriptorPortHamiltonian::new(n, n, 1, j, r, g, Box::new(storage.clone()))?;
        Ok(Self {
            phs,
            storage,
            profile,
            x: vec![0.0; n],
            card,
            rho: gas.density,
            moving_mass_kg,
            has_panel: panel.is_some(),
            time_s: 0.0,
        })
    }

    /// Joint electrical/mechanical step under held terminal voltage. Returns
    /// the native work/heat record and the same compact 1-m pressure observer
    /// as TsDriver. An invalid trial, unconverged solve or excursion refusal
    /// leaves the accepted driver state unchanged.
    pub fn step(
        &mut self,
        voltage_v: f64,
        duration_s: f64,
    ) -> Result<(StepRecord, f64), SpeakerError> {
        let target_s = self.time_s + duration_s;
        if !target_s.is_finite() || (duration_s > 0.0 && target_s <= self.time_s) {
            return Err(invalid(
                "magnetic speaker clock increment is unrepresentable",
            ));
        }
        let record = step_descriptor(&self.phs, &self.x, &[voltage_v], duration_s)?;
        self.profile.sample(record.x[1])?;
        if record.x[1].abs() > self.card.x_max_m {
            return Err(SpeakerError::ExcursionExceeded {
                excursion_m: record.x[1],
                x_max_m: self.card.x_max_m,
            });
        }
        let acceleration = if duration_s == 0.0 {
            0.0
        } else {
            (record.x[2] - self.x[2]) / (self.moving_mass_kg * duration_s)
        };
        let pressure = self.rho * self.card.sd_m2 * acceleration / (4.0 * core::f64::consts::PI);
        if !pressure.is_finite() {
            return Err(invalid("magnetic speaker pressure is unrepresentable"));
        }
        self.x.clone_from(&record.x);
        self.time_s = target_s;
        Ok((record, pressure))
    }

    // A thermal trial changes only the coil's dissipative stamp. It does not
    // reinitialize linkage, reinterpret stored energy or alter Bl/L/mechanics.
    // The accepted driver is immutable until the joint thermal candidate is
    // published. At most five states are copied for this reference owner.
    fn with_coil_resistance(&self, resistance_ohm: f64) -> Result<Self, SpeakerError> {
        if !resistance_ohm.is_finite() || resistance_ohm <= 0.0 {
            return Err(invalid("finite positive voice-coil resistance required"));
        }
        let n = self.phs.state_dim();
        let (j, old_r, g) = self.phs.structure();
        let mut r = old_r.to_vec();
        r[0] = resistance_ohm;
        let phs = DescriptorPortHamiltonian::new(
            n,
            n,
            1,
            j.to_vec(),
            r,
            g.to_vec(),
            Box::new(self.storage.clone()),
        )?;
        Ok(Self {
            phs,
            storage: self.storage.clone(),
            profile: self.profile,
            x: self.x.clone(),
            card: self.card,
            rho: self.rho,
            moving_mass_kg: self.moving_mass_kg,
            has_panel: self.has_panel,
            time_s: self.time_s,
        })
    }

    /// Accepted physical clock [s]. A refused step leaves this unchanged.
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Resistance used by the accepted electrical model [ohm]. The thermal
    /// adapter refreshes it from the actual material temperature after every
    /// jointly accepted step; standalone stepping holds it fixed.
    pub fn coil_resistance_ohm(&self) -> f64 {
        self.phs.structure().1[0]
    }

    /// Instantaneous physical voice-coil current [A].
    pub fn coil_current_a(&self) -> f64 {
        let local = self
            .profile
            .sample(self.x[1])
            .expect("accepted magnetic displacement");
        (self.x[0] - local.permanent_linkage_wb) / local.inductance_h
    }

    /// Instantaneous permanent-field plus reluctance force on the cone [N].
    pub fn magnetic_force_n(&self) -> f64 {
        let local = self
            .profile
            .sample(self.x[1])
            .expect("accepted magnetic displacement");
        let current = self.coil_current_a();
        current * local.force_factor_tm
            + 0.5 * current * current * local.inductance_gradient_h_per_m
    }

    /// Accepted cone displacement [m].
    pub fn excursion_m(&self) -> f64 {
        self.x[1]
    }

    /// Total magnetic, mechanical and cavity stored energy [J].
    pub fn stored_energy_j(&self) -> f64 {
        self.phs.hamiltonian(&self.x)
    }

    /// Whether the original cabinet-panel mode is composed in.
    pub const fn has_panel(&self) -> bool {
        self.has_panel
    }
}

#[cfg(test)]
mod tests;
