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

/// Bias-dependent small-signal analysis of the original reversible motor energy.
pub mod resonance {
    //! With i=(lambda-psi(x))/L(x), B=psi'(x), and C=B+i*L', the magnetic
    //! Hessian is [1/L, -C/L; -C/L, C^2/L-i*B'-(i^2/2)*L'']. These cross
    //! terms retain reciprocal back EMF, reluctance force and magnetic softening.
    //! Freezing only Bl/Le in a replacement TS card would lose the bias curvature.
    use super::{CoilStorage, MagneticTsDriver, Storage};
    use fs_phs::harmonic::equilibrium::{
        EquilibriumLinearization, EquilibriumPoint, LinearizationBudget, LinearizationError,
    };

    impl CoilStorage {
        fn magnetic_curvature(&self, state: &[f64]) -> Option<[f64; 3]> {
            if state.len()<3 || state.iter().any(|x|!x.is_finite()) {return None;}
            let x=state[1];
            let local=self.profile.sample(x).ok()?;
            let scale=self.profile.displacement_scale_m;
            let xi=x/scale;
            let current=(state[0]-local.permanent_linkage_wb)/local.inductance_h;
            let b_gradient=(self.profile.force_factor_tm[1]+2.0*xi*self.profile.force_factor_tm[2])/scale;
            let l_curvature=(2.0*self.profile.inductance_h[2]/scale)/scale;
            let coupled=local.force_factor_tm+current*local.inductance_gradient_h_per_m;
            let q00=1.0/local.inductance_h;
            let q01=-coupled/local.inductance_h;
            let q11=(coupled/local.inductance_h)*coupled-current*b_gradient
                -0.5*current*current*l_curvature;
            let block=[q00,q01,q11];
            if !current.is_finite() || block.iter().any(|x|!x.is_finite()) {None} else {Some(block)}
        }

        fn apply_curvature(&self, block:[f64;3], direction:&[f64], out:&mut[f64])->bool {
            let n=direction.len();
            if n<3 || out.len()!=n || n.checked_mul(n)!=Some(self.mechanical.matrix().len())
                || direction.iter().any(|v|!v.is_finite()) {return false;}
            self.mechanical.gradient(direction,out);
            out[0]+=block[0]*direction[0]+block[1]*direction[1];
            out[1]+=block[1]*direction[0]+block[2]*direction[1];
            out.iter().all(|x|x.is_finite())
        }
    }

    impl MagneticTsDriver {
        /// Read-only accepted [total linkage, cone displacement, cone momentum,
        /// optional panel displacement, optional panel momentum]. This is not
        /// [current, displacement, velocity]. No accepted state is changed.
        pub fn accepted_state(&self)->&[f64] {&self.x}

        /// Small-signal resonances and responses at a resolved original equilibrium.
        /// point.input is the held terminal voltage; state follows accepted_state().
        /// The source's actual resistance (including an admitted thermal refresh),
        /// complete mechanical storage, sealed cavity and panel coupling are retained.
        /// This does not solve for the bias or reset the accepted driver trajectory.
        ///
        /// Displacement must be strictly inside the source profile interval for a
        /// two-sided neighborhood. The original stationary equations must pass.
        /// Negative incremental stiffness is retained, not repaired to stability.
        /// Finite-amplitude harmonics, thermal perturbation dynamics, nonlinear
        /// suspension, saturation and moving acoustic geometry are not supplied.
        pub fn linearize_equilibrium(
            &self, point:EquilibriumPoint<'_>, budget:LinearizationBudget,
            mut cancelled:impl FnMut()->bool,
        )->Result<EquilibriumLinearization,LinearizationError> {
            if cancelled() {return Err(LinearizationError::Cancelled);}
            let limits=self.profile.displacement_limits_m;
            if point.state.len()!=self.phs.state_dim() || point.input.len()!=1
                || !point.state[1].is_finite() || point.state[1]<=limits[0] || point.state[1]>=limits[1] {
                return Err(LinearizationError::Invalid {
                    what:"complete magnetic-driver coordinates and an interior profile displacement required",
                });
            }
            let block=self.storage.magnetic_curvature(point.state).ok_or(LinearizationError::Invalid {
                what:"original magnetic-energy curvature is not representable",
            })?;
            // One constitutive evaluation for all columns; source and point are immutable.
            EquilibriumLinearization::from_descriptor(&self.phs,point,budget,
                &|_,direction,out|self.storage.apply_curvature(block,direction,out),&mut cancelled)
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use super::super::{PanelMode,TsCard,VoiceCoilMagnetics};
        use fs_math::c64::C64;
        use fs_phs::harmonic::{HarmonicBudget,HarmonicModel,PoleOptions};
        use std::cell::Cell;

        fn card()->TsCard {TsCard {re_ohm:0.8,le_h:0.6,bl:1.2,mms_kg:0.7,
            cms_m_per_n:0.9,rms_n_s_m:0.1,sd_m2:0.001,x_max_m:0.1}}
        fn profile(nonlinear:bool)->VoiceCoilMagnetics {
            VoiceCoilMagnetics::new([-0.05,0.05],0.01,
                if nonlinear {[1.2,-0.4,-0.2]} else {[1.2,0.,0.]},
                if nonlinear {[0.6,0.15,0.08]} else {[0.6,0.,0.]},).unwrap()
        }
        fn budget()->LinearizationBudget {LinearizationBudget {
            harmonic:HarmonicBudget {max_dimension:16,max_dense_bytes:2_000_000,
                max_dense_work:10_000_000,relative_residual_tolerance:1e-11},max_hessian_actions:16,
        }}
        fn close(a:C64,b:C64,tol:f64) {assert!((a-b).abs()<tol*(1.+b.abs()),"{a:?} vs {b:?}");}

        // Resolve a manufactured DC fixture from original static force and panel balance.
        fn bias(driver:&MagneticTsDriver,displacement:f64)->(Vec<f64>,f64,f64) {
            let n=driver.phs.state_dim();let q=driver.storage.mechanical.matrix();
            let mut state=vec![0.;n];state[1]=displacement;
            if driver.has_panel {state[3]=-q[3*n+1]*displacement/q[3*n+3];}
            let force=(0..n).map(|j|q[n+j]*state[j]).sum::<f64>();
            let local=driver.profile.sample(displacement).unwrap();
            let disc=(local.force_factor_tm.powi(2)+2.*local.inductance_gradient_h_per_m*force).sqrt();
            let current=2.*force/(local.force_factor_tm+disc);
            state[0]=local.permanent_linkage_wb+local.inductance_h*current;
            (state,driver.coil_resistance_ohm()*current,current)
        }
        fn snapshot(driver:&MagneticTsDriver,x:f64)->EquilibriumLinearization {
            let (state,voltage,_)=bias(driver,x);
            driver.linearize_equilibrium(EquilibriumPoint {state:&state,input:&[voltage],
                balance_scales:&vec![1.;state.len()],relative_tolerance:1e-10},budget(),||false).unwrap()
        }

        #[test]
        fn g1_original_energy_hessian_retains_reciprocal_and_curvature_terms() {
            let driver=MagneticTsDriver::new(card(),None,None,profile(true)).unwrap();
            for x in [-0.007,0.,0.009] {for current in [-0.2,0.,0.3] {
                let local=driver.profile.sample(x).unwrap();
                let state=[local.permanent_linkage_wb+local.inductance_h*current,x,0.1];
                let block=driver.storage.magnetic_curvature(&state).unwrap();
                for j in 0..3 {
                    let mut direction=[0.;3];direction[j]=1.;let mut exact=[0.;3];
                    assert!(driver.storage.apply_curvature(block,&direction,&mut exact));
                    let h=if j==1 {1e-8} else {1e-7};
                    let mut a=state;let mut b=state;a[j]+=h;b[j]-=h;
                    let mut ga=[0.;3];let mut gb=[0.;3];
                    driver.storage.gradient(&a,&mut ga);driver.storage.gradient(&b,&mut gb);
                    for i in 0..3 {
                        let difference=(ga[i]-gb[i])/(2.*h);
                        assert!((difference-exact[i]).abs()<2e-7*(1.+exact[i].abs()));
                    }
                }
            }}
        }

        #[test]
        fn g1_bias_response_matches_independent_current_coordinate_electromechanics() {
            let driver=MagneticTsDriver::new(card(),None,None,profile(true)).unwrap();
            let mut sampled=Vec::new();
            for x in [-0.005,0.,0.005] {
                let (_,_,current)=bias(&driver,x);let local=driver.profile.sample(x).unwrap();
                let scale=driver.profile.displacement_scale_m;let xi=x/scale;
                let bp=(driver.profile.force_factor_tm[1]+2.*xi*driver.profile.force_factor_tm[2])/scale;
                let lpp=2.*driver.profile.inductance_h[2]/(scale*scale);
                let coupled=local.force_factor_tm+current*local.inductance_gradient_h_per_m;
                let k=driver.storage.mechanical.matrix()[4]-current*bp-0.5*current*current*lpp;
                let mass=driver.moving_mass_kg;let damping=driver.phs.structure().1[8];
                let tangent=snapshot(&driver,x);
                assert_eq!(tangent.equilibrium_state().len(),3);
                assert!((tangent.equilibrium_output()[0]-current).abs()<1e-13);
                for w in [0.3,1.,3.,10.] {
                    let s=C64::new(0.,-w);let voltage=C64::new(0.2,-0.1);
                    let mechanics=(s*s).scale(mass)+s.scale(damping)+C64::from_re(k);
                    let z=s.scale(local.inductance_h)+C64::from_re(driver.coil_resistance_ohm());
                    let denominator=z*mechanics+s.scale(coupled*coupled);
                    let displacement=voltage.scale(coupled)/denominator;
                    let amps=voltage*mechanics/denominator;
                    let r=tangent.model().factor(w,||false).unwrap().response(&[voltage],||false).unwrap();
                    close(r.state[0],amps.scale(local.inductance_h)+displacement.scale(coupled),1e-11);
                    close(r.state[1],displacement,1e-11);
                    close(r.state[2],s.scale(mass)*displacement,1e-11);
                    close(r.output[0],amps,1e-11);
                }
                sampled.push(tangent.tangent_storage_matrix()[4]);
            }
            assert!((sampled[0]-sampled[1]).abs()>1e-3 && (sampled[2]-sampled[1]).abs()>1e-3);
        }

        #[test]
        fn g1_all_three_biased_poles_match_independent_current_coordinate_generator() {
            let driver=MagneticTsDriver::new(card(),None,None,profile(true)).unwrap();
            for x in [-0.005,0.,0.005] {
                let (_,_,current)=bias(&driver,x);let sample=driver.profile.sample(x).unwrap();
                let length=driver.profile.displacement_scale_m;let xi=x/length;
                let bp=(driver.profile.force_factor_tm[1]+2.*xi*driver.profile.force_factor_tm[2])/length;
                let lpp=2.*driver.profile.inductance_h[2]/(length*length);
                let coupled=sample.force_factor_tm+current*sample.inductance_gradient_h_per_m;
                let k=driver.storage.mechanical.matrix()[4]-current*bp-0.5*current*current*lpp;
                let mass=driver.moving_mass_kg;let loss=driver.phs.structure().1[8];
                let r=driver.coil_resistance_ohm();let l=sample.inductance_h;
                let a=[-r/l,0.,-coupled/l,0.,0.,1.,coupled/mass,-k/mass,-loss/mass].map(C64::from_re);
                let reference=fs_la::eigen_complex::eig(&a,3).unwrap();let local=snapshot(&driver,x);
                for s in reference {
                    let amps=((s*s).scale(mass)+s.scale(loss)+C64::from_re(k)).scale(1./coupled);
                    let seed=[amps.scale(l)+C64::from_re(coupled),C64::ONE,s.scale(mass)];
                    let pole=local.model().resonance_pole(s,&seed,&PoleOptions {
                        spectral_scale:5.,step_tolerance:1e-12,residual_tolerance:1e-12,max_iterations:32,
                        max_dimension:16,max_dense_bytes:2_000_000,max_work:10_000_000,
                    },||false).unwrap();
                    close(pole.laplace_pole_per_s,s,1e-10);
                    assert_eq!(pole.state.len(),3);
                }
            }
        }

        #[test]
        fn g3_original_cavity_and_panel_couplings_are_retained_at_zero_bias() {
            let panel=PanelMode {frequency_hz:2.,mass_kg:0.05,area_m2:0.001,damping_ratio:0.02};
            let nonlinear=MagneticTsDriver::new(card(),Some(0.02),Some(panel),profile(true)).unwrap();
            let constant=MagneticTsDriver::new(card(),Some(0.02),Some(panel),profile(false)).unwrap();
            let tangent=snapshot(&nonlinear,0.);
            let original=HarmonicModel::from_descriptor(&constant.phs,budget().harmonic,||false).unwrap();
            assert_eq!(tangent.tangent_storage_matrix().len(),25);
            assert!(tangent.tangent_storage_matrix()[1*5+3].abs()>0.);
            for w in [0.5,2.,6.,12.] {
                let a=tangent.model().factor(w,||false).unwrap().response(&[C64::ONE],||false).unwrap();
                let b=original.factor(w,||false).unwrap().response(&[C64::ONE],||false).unwrap();
                assert_eq!(a.state.len(),5);
                for (&x,&y) in a.state.iter().zip(&b.state) {close(x,y,1e-11);}
                close(a.output[0],b.output[0],1e-11);
            }
        }

        #[test]
        fn g0_g4_source_domain_equilibrium_and_late_cancellation_leave_driver_unchanged() {
            let driver=MagneticTsDriver::new(card(),None,None,profile(true)).unwrap();
            let before=driver.accepted_state().to_vec();let time=driver.time_s();
            let (state,voltage,_)=bias(&driver,0.005);let input=[voltage];let scales=[1.;3];
            let point=EquilibriumPoint {state:&state,input:&input,balance_scales:&scales,relative_tolerance:1e-10};
            let polls=Cell::new(0);
            let reference=driver.linearize_equilibrium(point,budget(),||{polls.set(polls.get()+1);false}).unwrap();
            let total=polls.get();
            for cut in [1,total/2,total] {
                let count=Cell::new(0);
                assert!(driver.linearize_equilibrium(point,budget(),||{count.set(count.get()+1);count.get()>=cut}).is_err());
            }
            let bad_input=[voltage+0.01];
            assert!(driver.linearize_equilibrium(EquilibriumPoint{input:&bad_input,..point},budget(),||false).is_err());
            for x in [-0.05,0.05,0.051] {
                let mut bad=state.clone();bad[1]=x;
                assert!(driver.linearize_equilibrium(EquilibriumPoint{state:&bad,..point},budget(),||false).is_err());
            }
            let replay=driver.linearize_equilibrium(point,budget(),||false).unwrap();
            assert_eq!(reference.tangent_storage_matrix(),replay.tangent_storage_matrix());
            assert_eq!(before.as_slice(),driver.accepted_state());assert_eq!(time,driver.time_s());
        }
    }
}
