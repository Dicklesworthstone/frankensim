//! Physical frame loads and observations on the shared modal/time/coupling owners.
//!
//! The original `Frame2d` owns the signed tangent, consistent reference mass,
//! condensation and eigenvectors. This adapter keeps their normalization:
//! `u = u_ref + Phi*q`, `Q = Phi^T*f - Gamma*a_base`. Thus nodal forces and
//! moments do the same work as their generalized forces. There is no second
//! eigensolver, time integrator, invented radiation transfer or damping fit.
//!
//! This is small vibration about a frozen, externally maintained reference
//! state, not nonlinear continuation of its material/geometry history. Loads
//! are increments about that reference; do not apply its gravity preload twice.
//! Omitted modes and supported inertia remain explicit in the mass inventory.
//! Outputs lie in the retained modal subspace: no residual static flexibility
//! or instantaneous load-dependent massless-coordinate feedthrough is inferred.

mod runtime;
pub use runtime::{FrameResponseLimits, FrameVibrationResponse, FrameVibrationRuntime};
mod material;
pub use material::FrameVibrationComponent;
pub use material::thermal;
pub use material::thermal::thermoelastic;

use fs_math::c64::C64;
use fs_solid::frame2d::{Frame2d, FrameMaterialDamping, FrameMode, FrameModeStatus};
use fs_solid::SolidError;

use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeError, ModalAcousticTimeModel,
};
use crate::render::schedule::force::coupled::ModalAttachment;

/// Admission before the original dense frame assembly/eigensolve allocates.
#[derive(Clone, Copy, Debug)]
pub struct FrameVibrationBudget {
    /// Full physical coordinates, including fixed and massless coordinates.
    pub max_frame_dofs: usize,
    /// Full `ndof*ndof` entries, not merely the retained modal output.
    /// The native frame/eigensolver has its own additional admission limits.
    pub max_dense_entries: usize,
}

/// A refused frame-to-modal binding or physical projection.
#[derive(Debug)]
pub enum FrameVibrationError {
    /// Invalid dimensions, physical domain, budget or unrepresentable result.
    Invalid(&'static str),
    /// Original structural assembly or eigensolver refusal.
    Solid(SolidError),
    /// Original shared time-domain owner refusal.
    Modal(ModalAcousticTimeError),
    /// Original simultaneous modal connection/propagation owner refusal.
    Coupling(crate::render::schedule::force::coupled::ModalCouplingError),
    /// No unstable mode may be hidden by retaining a smaller positive subset.
    UnstableMode { mode: usize, eigenvalue_per_s2: f64 },
    /// A near-zero mode cannot be silently converted into a positive oscillator.
    UnresolvedMode { mode: usize, eigenvalue_per_s2: f64 },
    /// Cancellation before publication; output scratch may have been overwritten.
    Cancelled,
}
impl core::fmt::Display for FrameVibrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frame vibration: {self:?}")
    }
}
impl std::error::Error for FrameVibrationError {}
impl From<SolidError> for FrameVibrationError {
    fn from(error: SolidError) -> Self { Self::Solid(error) }
}
impl From<ModalAcousticTimeError> for FrameVibrationError {
    fn from(error: ModalAcousticTimeError) -> Self { Self::Modal(error) }
}
impl From<crate::render::schedule::force::coupled::ModalCouplingError> for FrameVibrationError {
    fn from(error: crate::render::schedule::force::coupled::ModalCouplingError) -> Self {
        match error {
            crate::render::schedule::force::coupled::ModalCouplingError::Cancelled => Self::Cancelled,
            other => Self::Coupling(other),
        }
    }
}

/// An existing modal mechanical port plus its frozen physical reference datum.
#[derive(Clone, Debug)]
pub struct FrameTranslationAttachment {
    /// The original reciprocal force/displacement map, in incremental coordinates.
    pub port: ModalAttachment,
    /// Add to the port displacement to recover total reference-relative position.
    /// Account for this datum when specifying spring rest extensions/contact gaps.
    pub reference_displacement_m: f64,
}

/// Immutable source-derived basis, shared by forcing, observations and contacts.
#[derive(Clone, Debug)]
pub struct FrameVibrationBasis {
    modes: Vec<FrameMode>,
    reference_displacement: Vec<f64>,
    available_modes: usize,
    influence_mass_kg: [f64; 2],
    effective_mass_kg: [f64; 2],
    material_damping: Option<FrameMaterialDamping>,
}

impl FrameVibrationBasis {
    /// Extract the complete signed spectrum, refuse any unstable/unresolved
    /// mode, then explicitly retain the lowest `retained_modes` oscillators.
    /// No Nyquist-based or stability-based silent truncation occurs.
    ///
    /// Cancellation is checked around the bounded native eigensolve and per
    /// mode; this adapter cannot preempt that synchronous eigensolve itself.
    /// Existing dynamic velocity is NOT imported: the returned time model is
    /// initially at rest about `frame.displacements()`.
    pub fn from_frame(
        frame: &Frame2d,
        retained_modes: usize,
        budget: FrameVibrationBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FrameVibrationError> {
        poll(&mut cancelled)?;
        let n = frame.ndof();
        if n == 0 || n > budget.max_frame_dofs || retained_modes == 0
            || n.checked_mul(n).is_none_or(|entries| entries > budget.max_dense_entries)
        {
            return Err(FrameVibrationError::Invalid("full frame/dense preparation budget or retained mode count"));
        }
        let mut analysis = frame.modal_analysis()?;
        poll(&mut cancelled)?;
        if retained_modes > analysis.modes.len() {
            return Err(FrameVibrationError::Invalid("retained count exceeds the complete frame spectrum"));
        }
        for (index, mode) in analysis.modes.iter().enumerate() {
            poll(&mut cancelled)?;
            match mode.status() {
                FrameModeStatus::Oscillatory => {},
                FrameModeStatus::Unstable => return Err(FrameVibrationError::UnstableMode {
                    mode: index, eigenvalue_per_s2: mode.eigenvalue_per_s2,
                }),
                FrameModeStatus::Unresolved => return Err(FrameVibrationError::UnresolvedMode {
                    mode: index, eigenvalue_per_s2: mode.eigenvalue_per_s2,
                }),
            }
        }
        let available_modes = analysis.modes.len();
        analysis.modes.truncate(retained_modes);
        let effective_mass_kg = [
            checked_sum(analysis.modes.iter().map(|m| m.effective_mass_kg[0]))?,
            checked_sum(analysis.modes.iter().map(|m| m.effective_mass_kg[1]))?,
        ];
        let reference_displacement = frame.displacements().to_vec();
        if reference_displacement.iter().any(|x| !x.is_finite()) {
            return Err(FrameVibrationError::Invalid("nonfinite frozen frame displacement"));
        }
        // The same immutable source operator serves native Newmark and every
        // reduced receiver. It retains the complete spectrum's sampling gate,
        // so truncation cannot hide an unresolved material-damped frame mode.
        let material_damping = frame.material_damping_operator()?;
        poll(&mut cancelled)?;
        Ok(Self {
            modes: analysis.modes, reference_displacement, available_modes,
            influence_mass_kg: analysis.base_influence_mass_kg, effective_mass_kg,
            material_damping,
        })
    }

    /// Retained source eigenpairs, with their residual and resolution diagnostics.
    pub fn modes(&self) -> &[FrameMode] { &self.modes }
    /// Full physical coordinate count, node-major `[ux, uy, theta]`.
    pub fn ndof(&self) -> usize { self.reference_displacement.len() }
    /// Complete stable source spectrum count before explicit truncation.
    pub fn available_mode_count(&self) -> usize { self.available_modes }
    /// Frozen accepted displacement; dynamic outputs are increments about it.
    pub fn reference_displacement(&self) -> &[f64] { &self.reference_displacement }
    /// Full horizontal/vertical reference influence inventory, including support mass.
    pub fn base_influence_mass_kg(&self) -> [f64; 2] { self.influence_mass_kg }
    /// Sum of retained squared participation factors, not the complete inventory.
    pub fn retained_effective_mass_kg(&self) -> [f64; 2] { self.effective_mass_kg }
    /// Original physical viscosity operator and all retained member sources.
    pub fn material_damping(&self) -> Option<&FrameMaterialDamping> {
        self.material_damping.as_ref()
    }

    /// Construct the ORIGINAL exact held-force modal owner, with explicitly
    /// supplied per-mode viscous damping ratios. No acoustic pressure transfer
    /// is guessed from a planar structural model. This component can also be
    /// passed directly to `CoupledModalSystem` and its existing contact wrappers.
    /// Material-viscous frames require [`Self::coupled_component`], which carries
    /// their complete nonproportional damping into that SAME coupled owner.
    /// This independent-mode path refuses rather than dropping those losses.
    pub fn time_model(
        &self,
        sample_rate_hz: u32,
        damping_ratios: &[f64],
        budget: ModalAcousticTimeBudget,
    ) -> Result<ModalAcousticTimeModel, FrameVibrationError> {
        if self.material_damping.is_some() {
            return Err(FrameVibrationError::Invalid(
                "material frame viscosity requires the complete coupled_component binding",
            ));
        }
        self.oscillator_model(sample_rate_hz, damping_ratios, budget)
    }

    fn oscillator_model(
        &self,
        sample_rate_hz: u32,
        damping_ratios: &[f64],
        budget: ModalAcousticTimeBudget,
    ) -> Result<ModalAcousticTimeModel, FrameVibrationError> {
        if damping_ratios.len() != self.modes.len() {
            return Err(FrameVibrationError::Invalid("one explicit damping ratio per retained mode is required"));
        }
        let modes = self.modes.iter().zip(damping_ratios).map(|(mode, &damping_ratio)| {
            ModalAcousticMode {
                angular_frequency_rad_s: mode.eigenvalue_per_s2.sqrt(),
                damping_ratio,
                pressure_per_modal_velocity: C64::ZERO,
            }
        }).collect();
        Ok(ModalAcousticTimeModel::try_new(sample_rate_hz, modes, budget)?)
    }

    /// Project incremental nodal forces [N] and moments [N m], together with
    /// simultaneous horizontal/vertical uniform base acceleration [m/s^2].
    /// `Q = Phi^T*f - Gamma_x*a_x - Gamma_y*a_y` in N/sqrt(kg).
    /// Relative-frame inertial work is not absolute actuator work.
    ///
    /// No allocation; `generalized` is scratch and may be partly overwritten
    /// on overflow/cancellation. No accepted physical state is changed.
    pub fn project_force_into(
        &self,
        nodal_force_and_moment: &[f64],
        base_acceleration_m_s2: [f64; 2],
        generalized: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), FrameVibrationError> {
        poll(&mut cancelled)?;
        if nodal_force_and_moment.len() != self.ndof() || generalized.len() != self.modes.len()
            || nodal_force_and_moment.iter().chain(base_acceleration_m_s2.iter()).any(|x| !x.is_finite())
        {
            return Err(FrameVibrationError::Invalid("finite physical loads, base acceleration and matching modal output required"));
        }
        for (mode, force) in self.modes.iter().zip(generalized) {
            poll(&mut cancelled)?;
            let applied = checked_sum(mode.displacement.iter().zip(nodal_force_and_moment).map(|(p, f)| p*f))?;
            *force = checked_sum([
                applied,
                -mode.base_participation[0]*base_acceleration_m_s2[0],
                -mode.base_participation[1]*base_acceleration_m_s2[1],
            ].into_iter())?;
        }
        poll(&mut cancelled)
    }

    /// Reconstruct total reference-relative displacement and relative velocity.
    /// Translational entries use m and m/s; angular entries use rad and rad/s.
    /// The base's own translation/velocity is deliberately not added.
    /// Output buffers are scratch on a refusal. No allocation occurs.
    pub fn reconstruct_into(
        &self,
        states: &[ModalAcousticState],
        displacement: &mut [f64],
        velocity: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), FrameVibrationError> {
        poll(&mut cancelled)?;
        if states.len() != self.modes.len() || displacement.len() != self.ndof() || velocity.len() != self.ndof()
            || states.iter().any(|s| !s.displacement_m_sqrt_kg.is_finite() || !s.velocity_m_sqrt_kg_per_s.is_finite())
        {
            return Err(FrameVibrationError::Invalid("finite modal state and complete physical output buffers required"));
        }
        for i in 0..self.ndof() {
            if i % 128 == 0 { poll(&mut cancelled)?; }
            displacement[i] = checked_sum(core::iter::once(self.reference_displacement[i]).chain(
                self.modes.iter().zip(states).map(|(mode, state)| mode.displacement[i]*state.displacement_m_sqrt_kg)))?;
            velocity[i] = checked_sum(self.modes.iter().zip(states)
                .map(|(mode, state)| mode.displacement[i]*state.velocity_m_sqrt_kg_per_s))?;
        }
        poll(&mut cancelled)
    }

    /// Mechanical attachment at an original node along an explicitly unit
    /// translation direction. The SAME shape maps contact force into modes and
    /// modal motion back to the attachment, preserving virtual work.
    /// These are incremental ports; the separate reference datum is not lost.
    pub fn translation_attachment(
        &self,
        component: usize,
        node: usize,
        direction: [f64; 2],
    ) -> Result<FrameTranslationAttachment, FrameVibrationError> {
        if node >= self.ndof()/3 || direction.iter().any(|x| !x.is_finite())
            || (direction[0].hypot(direction[1]) - 1.0).abs() > 1e-12
        {
            return Err(FrameVibrationError::Invalid("existing frame node and finite unit attachment direction required"));
        }
        let coordinate = |shape: &[f64]| checked_sum([
            direction[0]*shape[3*node], direction[1]*shape[3*node+1],
        ].into_iter());
        let shapes = self.modes.iter().map(|mode| coordinate(&mode.displacement))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(FrameTranslationAttachment {
            port: ModalAttachment { component, shapes },
            reference_displacement_m: coordinate(&self.reference_displacement)?,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FrameVibrationError> {
    if cancelled() { Err(FrameVibrationError::Cancelled) } else { Ok(()) }
}

/// Compensated deterministic reduction for physical power-conjugate maps.
fn checked_sum(values: impl Iterator<Item = f64>) -> Result<f64, FrameVibrationError> {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        if !value.is_finite() { return Err(FrameVibrationError::Invalid("nonfinite frame projection term")); }
        let next = sum + value;
        correction += if sum.abs() >= value.abs() { (sum-next)+value } else { (value-next)+sum };
        sum = next;
    }
    let result = sum + correction;
    if result.is_finite() { Ok(result) } else { Err(FrameVibrationError::Invalid("frame projection overflow")) }
}

#[cfg(test)]
mod runtime_tests {
    use crate::frame_vibration::{FrameResponseLimits, FrameVibrationBasis,
        FrameVibrationBudget, FrameVibrationError, FrameVibrationRuntime};
    use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
    use fs_solid::frame2d::{Frame2d, Geometry, SectionModel};

    fn limits() -> FrameResponseLimits {
        FrameResponseLimits {
            maximum_incremental_translation_m: 1.0, maximum_incremental_rotation_rad: 2.0,
            maximum_translation_speed_m_s: 10.0, maximum_angular_speed_rad_s: 20.0,
            maximum_absolute_translation_acceleration_m_s2: 100.0,
            maximum_angular_acceleration_rad_s2: 200.0,
        }
    }
    fn basis(preload: f64) -> FrameVibrationBasis {
        let mut frame = Frame2d::new(Geometry::Linear);
        frame.add_node(0.0, 0.0); frame.add_node(1.0, 0.0);
        frame.fix(0, [true; 3]); frame.fix(1, [true, false, false]);
        frame.add_mass(1, 2.0);
        frame.add_element(0, 1, 3, 0.0, &|| SectionModel::Elastic { ea: 200.0, ei: 8.0 }).unwrap();
        if preload != 0.0 { frame.static_load(&[0.0,0.0,0.0,0.0,preload,0.0], 1).unwrap(); }
        FrameVibrationBasis::from_frame(&frame, 1,
            FrameVibrationBudget { max_frame_dofs: 6, max_dense_entries: 36 }, || false).unwrap()
    }
    fn runtime(zeta: f64, initial: ModalAcousticState, limit: FrameResponseLimits) -> FrameVibrationRuntime {
        FrameVibrationRuntime::new(basis(0.0), 1000, &[zeta],
            ModalAcousticTimeBudget::audible_reference(), limit, &[initial], || false).unwrap()
    }
    fn close(a:f64,b:f64) { assert!((a-b).abs() < 2e-9*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}"); }

    #[test]
    fn held_nodal_load_reaches_analytic_motion_with_closed_work() {
        let mut sim = runtime(0.0, ModalAcousticState::default(), limits());
        let (mut work, mut loss) = (0.0, 0.0);
        for _ in 0..125 {
            let r = sim.step(&[0.0,0.0,0.0,0.0,4.0,0.0], [0.0; 2], || false).unwrap();
            work += r.input_work_j; loss += r.viscous_dissipation_j;
            close(r.nodal_work_j, r.input_work_j);
            assert_eq!(r.base_inertial_work_j, 0.0);
            assert!(r.work_projection_residual_j.abs() <= r.work_projection_tolerance_j);
        }
        let r = sim.response();
        let omega = 12.0_f64.sqrt();
        assert_eq!(r.sample, 125); assert_eq!(r.time_s, 0.125);
        close(r.displacement[4], 4.0/24.0*(1.0-(omega*0.125).cos()));
        close(r.relative_velocity[4], 4.0/24.0*omega*(omega*0.125).sin());
        close(r.relative_acceleration[4], 2.0*(omega*0.125).cos());
        close(r.retained_energy_j+loss, work);
        close(work, 4.0*r.displacement[4]);
    }

    #[test]
    fn accelerometer_output_includes_base_and_work_remains_relative() {
        let zeta = 0.03;
        let mut sim = runtime(zeta, ModalAcousticState::default(), limits());
        for _ in 0..40 {
            let r = sim.step(&[0.0,0.0,0.0,0.0,3.0,0.0], [0.2,0.4], || false).unwrap();
            close(r.absolute_acceleration[0], 0.2);
            close(r.absolute_acceleration[1], 0.4);
            close(r.relative_acceleration[0], 0.0);
            close(r.relative_acceleration[1], 0.0);
            close(r.absolute_acceleration[4]-r.relative_acceleration[4], 0.4);
            close(2.0*r.absolute_acceleration[4] + 4.0*zeta*12.0_f64.sqrt()*r.relative_velocity[4]
                +24.0*r.incremental_displacement[4], 3.0);
            close(r.input_work_j, r.nodal_work_j+r.base_inertial_work_j);
            close(r.base_inertial_work_j, -0.8/3.0*r.nodal_work_j);
        }
    }

    #[test]
    fn initial_pluck_ringdown_reuses_shared_damping_and_storage() {
        let initial = ModalAcousticState { displacement_m_sqrt_kg: 0.2, velocity_m_sqrt_kg_per_s: -0.1 };
        let mut sim = runtime(0.1, initial, limits());
        let energy0 = sim.response().retained_energy_j;
        close(energy0, 0.245);
        let (mut loss, mut previous) = (0.0, energy0);
        for _ in 0..500 {
            let r = sim.step(&[0.0; 6], [0.0; 2], || false).unwrap();
            assert!(r.retained_energy_j <= previous+1e-13);
            assert!(r.viscous_dissipation_j >= -r.dissipation_roundoff_tolerance_j);
            assert_eq!(r.input_work_j, 0.0);
            loss += r.viscous_dissipation_j;
            previous = r.retained_energy_j;
        }
        close(previous+loss, energy0);
        assert!(loss > 0.01);
    }

    #[test]
    fn late_physical_limit_preserves_modal_state_response_and_clock() {
        let mut bound = limits(); bound.maximum_incremental_translation_m = 1e-12;
        let mut sim = runtime(0.0, ModalAcousticState::default(), bound);
        let before = sim.response().clone();
        let states = sim.states().to_vec();
        assert!(sim.step(&[0.0,0.0,0.0,0.0,4.0,0.0], [0.0; 2], || false).is_err());
        assert_eq!(sim.response(), &before); assert_eq!(sim.states(), states.as_slice());
        let resumed = sim.step(&[0.0; 6], [0.0; 2], || false).unwrap();
        assert_eq!(resumed.sample, 1); assert_eq!(resumed.retained_energy_j, 0.0);
    }

    #[test]
    fn final_cancellation_rolls_back_and_retry_matches_uninterrupted_history() {
        let load = [0.0,0.0,0.0,0.0,4.0,0.0];
        let mut reference = runtime(0.02, ModalAcousticState::default(), limits());
        let mut checks = 0;
        reference.step(&load, [0.1,0.2], || { checks += 1; false }).unwrap();
        let mut sim = runtime(0.02, ModalAcousticState::default(), limits());
        let before = sim.response().clone();
        let states = sim.states().to_vec();
        let mut seen = 0;
        assert!(matches!(sim.step(&load, [0.1,0.2], || { seen += 1; seen == checks }),
            Err(FrameVibrationError::Cancelled)));
        assert_eq!(seen, checks);
        assert_eq!(sim.response(), &before); assert_eq!(sim.states(), states.as_slice());
        sim.step(&load, [0.1,0.2], || false).unwrap();
        assert_eq!(sim.response(), reference.response()); assert_eq!(sim.states(), reference.states());
    }

    #[test]
    fn preload_is_maintained_without_double_applying_its_force_or_energy() {
        let mut sim = FrameVibrationRuntime::new(basis(4.0), 1000, &[0.02],
            ModalAcousticTimeBudget::audible_reference(), limits(),
            &[ModalAcousticState::default()], || false).unwrap();
        for _ in 0..8 {
            let r = sim.step(&[0.0; 6], [0.0; 2], || false).unwrap();
            close(r.displacement[4], 4.0/24.0);
            assert_eq!(r.incremental_displacement, vec![0.0; 6]);
            assert_eq!(r.retained_energy_j, 0.0); assert_eq!(r.input_work_j, 0.0);
        }
    }

    #[test]
    fn invalid_initial_and_trial_data_cannot_publish_a_response() {
        let mut bad = limits(); bad.maximum_angular_speed_rad_s = f64::NAN;
        assert!(FrameVibrationRuntime::new(basis(0.0), 1000, &[0.0],
            ModalAcousticTimeBudget::audible_reference(), bad,
            &[ModalAcousticState::default()], || false).is_err());
        let mut sim = runtime(0.0, ModalAcousticState::default(), limits());
        let initial = sim.response().clone();
        assert!(sim.step(&[0.0; 5], [0.0; 2], || false).is_err());
        assert!(sim.step(&[0.0; 6], [f64::NAN,0.0], || false).is_err());
        assert!(sim.step(&[0.0; 6], [101.0,0.0], || false).is_err());
        assert_eq!(sim.response(), &initial);
    }

    #[test]
    fn physical_frames_exchange_vibration_through_the_original_reciprocal_network() {
        use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalConnection, ModalCouplingConfig};
        use fs_exec::CancelGate;
        let left = basis(0.0); let right = basis(0.0);
        let mut left_model = left.time_model(100, &[0.0], ModalAcousticTimeBudget::audible_reference()).unwrap();
        let right_model = right.time_model(100, &[0.0], ModalAcousticTimeBudget::audible_reference()).unwrap();
        left_model.restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.02/left.modes()[0].displacement[4],
            velocity_m_sqrt_kg_per_s: 0.0,
        }]).unwrap();
        let link = ModalConnection {
            left: left.translation_attachment(0, 1, [0.0,1.0]).unwrap().port,
            right: right.translation_attachment(1, 1, [0.0,1.0]).unwrap().port,
            stiffness_n_m: 6.0, damping_n_s_m: 0.0, rest_extension_m: 0.0,
        };
        let mut system = CoupledModalSystem::new(vec![left_model,right_model], vec![link], ModalCouplingConfig {
            max_modes: 2, max_connections: 1, max_setup_terms: 8192,
            nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
            maximum_abs_pressure_pa: 100.0, maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-10, energy_absolute_tolerance_j: 1e-11,
            energy_relative_tolerance: 1e-9,
        }, &CancelGate::new()).unwrap();
        close(system.total_energy_j().unwrap(), 0.006);
        for _ in 0..20 { system.step(&[0.0; 2]).unwrap(); }
        let (mut ul, mut ur, mut vl, mut vr) = ([0.0; 6],[0.0; 6],[0.0; 6],[0.0; 6]);
        left.reconstruct_into(system.components()[0].states(), &mut ul, &mut vl, || false).unwrap();
        right.reconstruct_into(system.components()[1].states(), &mut ur, &mut vr, || false).unwrap();
        let (symmetric, antisymmetric) = ((12.0_f64.sqrt()*0.2).cos(), (18.0_f64.sqrt()*0.2).cos());
        close(ul[4], 0.01*(symmetric+antisymmetric));
        close(ur[4], 0.01*(symmetric-antisymmetric));
        close(system.total_energy_j().unwrap(), 0.006);
        assert!(ur[4] > 0.001);
    }

    #[test]
    fn rigid_motion_and_compressive_instability_are_not_hidden_by_mode_selection() {
        use fs_solid::frame2d::FrameMassFormulation;
        let mut free = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
        free.add_node(0.0,0.0); free.add_node(1.0,0.0);
        free.add_element(0,1,3,2.0,&|| SectionModel::Elastic { ea: 200.0, ei: 8.0 }).unwrap();
        let admission = FrameVibrationBudget { max_frame_dofs: 6, max_dense_entries: 36 };
        assert!(matches!(FrameVibrationBasis::from_frame(&free,1,admission,||false),
            Err(FrameVibrationError::UnresolvedMode { .. })));
        let mut compressed = Frame2d::new(Geometry::Corotational);
        compressed.add_node(0.0,0.0); compressed.add_node(1.0,0.0);
        compressed.fix(0,[true;3]); compressed.add_mass(1,2.0);
        compressed.add_element(0,1,3,0.0,&|| SectionModel::Elastic { ea: 2.0e6, ei: 8.0 }).unwrap();
        compressed.static_load(&[0.0,0.0,0.0,-40.0,0.0,0.0],4).unwrap();
        assert!(matches!(FrameVibrationBasis::from_frame(&compressed,1,admission,||false),
            Err(FrameVibrationError::UnstableMode { .. })));
    }
}
