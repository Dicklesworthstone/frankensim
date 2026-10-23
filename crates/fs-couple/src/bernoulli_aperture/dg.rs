//! A moving aperture coupled to the original exact DG propagation owner.
//!
//! The DG exponential exposes mean pressure as `history + impedance * flow`.
//! The existing strict nonlinear aperture solver already accepts exactly this
//! law through `incoming = history / 2`. Its ExactZoh mechanical block supplies
//! mean swept flow from the complete opening displacement. The resulting mean
//! flow is held at the DG port; the original propagator integrates its pressure
//! and work. No fitted reflection, lagged pressure, second pressure solver, or
//! split wall-memory update enters the coupling.
//!
//! This is a bounded offline adapter. Candidate stepping allocates through the
//! original owners. Inlet pressure is an internal observation, not radiated
//! microphone pressure. Nonlinear held-force/time resolution and DG spatial
//! convergence remain separate from the linear exponential and Bessel-band
//! admissions. The dense DG loss ledger still infers loss from energy deficit.

use super::dynamic::{ApertureDrive, ApertureFrame, DynamicAperture};
use super::paired::BladePropagator;
use super::tube::TubeDrive;
use crate::acoustic_realize::AcousticRealizeError;
use fs_duct::dg::{DgBoundary, DgDuct, DgPropagator, DgState, DgStepOutput, PortId, PortInput};

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// Cold, source-bound one-port load. Read its mean impedance before constructing
/// the aperture, then move both into [`ApertureDg`]. Every matrix and wall state
/// remains owned by the original DG implementation.
#[derive(Debug, Clone)]
pub struct ApertureDgLoad {
    duct: DgDuct,
    propagator: DgPropagator,
    state: DgState,
    inlet: PortId,
    impedance_pa_s_m3: f64,
}

impl ApertureDgLoad {
    /// Prepare an initially quiescent duct with one inlet port and a closed,
    /// ideal-open, or matched outlet. `maximum_dofs` bounds the complete acoustic
    /// plus wall-memory state before constructing its dense exponential.
    ///
    /// # Errors
    /// Missing inlet, another unknown port, memory-budget excess, unresolved
    /// exponential, wall band beyond mechanical Nyquist, or nonpositive mean
    /// impedance. A zero lossless return-period response is refused, not repaired.
    pub fn new(
        duct: DgDuct,
        time_step_s: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
    ) -> Result<Self, AcousticRealizeError> {
        let DgBoundary::Port(inlet) = duct.inlet() else {
            return Err(invalid("DG aperture load requires one declared inlet port"));
        };
        if matches!(duct.outlet(), DgBoundary::Port(_)) {
            return Err(invalid(
                "DG aperture load does not independently solve a second unknown port",
            ));
        }
        if duct.total_dofs() > maximum_dofs {
            return Err(invalid(
                "DG aperture acoustic plus memory states exceed the declared DOF budget",
            ));
        }
        if let Some(model) = duct.viscothermal_model() {
            if model.band_rad_s()[1] * time_step_s > core::f64::consts::PI {
                return Err(invalid(
                    "DG aperture viscothermal band exceeds the mechanical Nyquist frequency",
                ));
            }
        }
        let propagator = duct
            .propagator(time_step_s, exponential_tolerance)
            .map_err(AcousticRealizeError::Duct)?;
        let state = duct.zero_state();
        let impedance_pa_s_m3 = propagator
            .port_response(&state, inlet)
            .map_err(AcousticRealizeError::Duct)?
            .impedance_pa_s_m3;
        if impedance_pa_s_m3 <= 0.0 {
            return Err(invalid(
                "DG aperture needs strictly positive exact mean-port impedance",
            ));
        }
        Ok(Self {
            duct,
            propagator,
            state,
            inlet,
            impedance_pa_s_m3,
        })
    }

    /// Exact-ZOH mean impedance to bind into the original aperture spec.
    #[must_use]
    pub const fn impedance_pa_s_m3(&self) -> f64 {
        self.impedance_pa_s_m3
    }

    /// Physical duct and gas/pole source, unchanged by the adapter.
    #[must_use]
    pub const fn duct(&self) -> &DgDuct {
        &self.duct
    }

    /// Original propagator, including its admitted step and certificate.
    #[must_use]
    pub const fn propagator(&self) -> &DgPropagator {
        &self.propagator
    }
}

/// One completely accepted physical step, including the original DG ledger.
#[derive(Debug, Clone, PartialEq)]
pub struct DgApertureFrame {
    /// Original nonlinear aperture observation and accepted mechanical state.
    pub aperture: ApertureFrame,
    /// Original exact propagation observation; its port pressure is a step mean.
    pub duct: DgStepOutput,
    /// Mechanical, contact, material, acoustic, and wall-memory storage [J].
    pub stored_energy_j: f64,
    /// Change in the complete storage [J].
    pub storage_change_j: f64,
    /// Aperture loss plus inferred duct loss [J]. This is NOT an independent
    /// wall/terminal heat integral when `duct.ledger.balance_checked` is false.
    pub inferred_total_dissipation_j: f64,
    /// Upstream pressure work on aperture flow, excluding body-supplied flow [J].
    pub upstream_work_j: f64,
    /// Inlet pressure work on independently supplied body flow [J].
    pub body_work_j: f64,
    /// DG integrated work minus aperture mean pressure times total bore volume
    /// [J]. The internal port cancels before either candidate is published.
    pub port_work_residual_j: f64,
}

impl DgApertureFrame {
    /// Diagnostic whole-system balance [J]. Inferred DG loss does not constitute
    /// a new independent energy certificate.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.storage_change_j + self.inferred_total_dissipation_j
            - self.upstream_work_j
            - self.body_work_j
            - self.aperture.mechanical_work_j
    }
}

/// One clock and one accepted state for the existing nonlinear aperture and DG
/// duct. Failed input, nonlinear solve, propagation or port compatibility leaves
/// both participants unchanged; source material/contact histories stay retained.
pub struct ApertureDg {
    aperture: DynamicAperture,
    load: ApertureDgLoad,
}

impl ApertureDg {
    /// Bind an unadvanced ExactZoh aperture to the prepared source load.
    ///
    /// # Errors
    /// A previously advanced aperture, a different gas density/clock/impedance,
    /// or a midpoint mechanical block. No participant is silently retuned.
    pub fn new(
        aperture: DynamicAperture,
        load: ApertureDgLoad,
    ) -> Result<Self, AcousticRealizeError> {
        let spec = aperture.spec();
        if aperture.accepted_steps() != 0
            || load.state.step_count() != 0
            || aperture.blade_propagator() != BladePropagator::ExactZoh
            || spec.time_step_s.to_bits() != load.propagator.dt().to_bits()
            || spec.density_kg_m3.to_bits() != load.duct.inlet_gas().density.to_bits()
            || spec.impedance_pa_s_m3.to_bits() != load.impedance_pa_s_m3.to_bits()
        {
            return Err(invalid(
                "DG coupling requires a new ExactZoh aperture bound to its mean impedance, gas density and clock",
            ));
        }
        Ok(Self { aperture, load })
    }

    /// Original retained aperture, including plate, closure, and material arms.
    #[must_use]
    pub const fn aperture(&self) -> &DynamicAperture {
        &self.aperture
    }

    /// Immutable cold source and original exponential metadata.
    #[must_use]
    pub const fn load(&self) -> &ApertureDgLoad {
        &self.load
    }

    /// Accepted acoustic and memory state; its physical nodal arrays are current.
    #[must_use]
    pub const fn duct_state(&self) -> &DgState {
        &self.load.state
    }

    /// Complete accepted stored energy [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.aperture.stored_energy_j() + self.load.state.stored_energy_j()
    }

    /// Increase the original aperture horizon without resetting either owner.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        self.aperture.extend_step_budget(total)
    }

    /// Advance with upstream pressure and an optional separately credited body
    /// flow. Both are held for this mechanical sample.
    pub fn step(&mut self, drive: TubeDrive) -> Result<DgApertureFrame, AcousticRealizeError> {
        self.step_with_force(drive, 0.0)
    }

    /// Solve the original aperture against the DG mean-port law, then validate
    /// candidate exact propagation and internal work before publishing either.
    ///
    /// # Errors
    /// Original strict nonlinear/admission errors, nonfinite observation, or
    /// inconsistent mean pressure. All failures preserve both accepted states.
    pub fn step_with_force(
        &mut self,
        drive: TubeDrive,
        force_n: f64,
    ) -> Result<DgApertureFrame, AcousticRealizeError> {
        let relation = self
            .load
            .propagator
            .port_response(&self.load.state, self.load.inlet)
            .map_err(AcousticRealizeError::Duct)?;
        let trial = self.aperture.preview_step_with_force(
            ApertureDrive {
                upstream_pressure_pa: drive.upstream_pressure_pa,
                incoming_pressure_pa: relation.history_pressure_pa / 2.0,
                body_flow_m3_s: drive.body_flow_m3_s,
            },
            force_n,
        )?;
        let aperture = trial.frame;
        let mut candidate = self.load.state.clone();
        let duct = self
            .load
            .propagator
            .step(
                &mut candidate,
                &[PortInput {
                    id: self.load.inlet,
                    effort: aperture.bore_flow_m3_s,
                }],
            )
            .map_err(AcousticRealizeError::Duct)?;
        let port = duct.port_outputs[0];
        let scale = relation.history_pressure_pa.abs()
            + (relation.impedance_pa_s_m3 * aperture.bore_flow_m3_s).abs()
            + aperture.bore_pressure_pa.abs();
        if !port.conjugate_rate.is_finite()
            || (port.conjugate_rate - aperture.bore_pressure_pa).abs() > 1e-10 * scale.max(1e-300)
        {
            return Err(invalid(
                "DG and aperture disagree on same-step mean port pressure",
            ));
        }
        let dt = self.load.propagator.dt();
        let duct_change = candidate.stored_energy_j() - self.load.state.stored_energy_j();
        let frame = DgApertureFrame {
            aperture,
            stored_energy_j: aperture.stored_energy_j + candidate.stored_energy_j(),
            storage_change_j: aperture.storage_change_j + duct_change,
            inferred_total_dissipation_j: aperture.dissipated_energy_j
                + duct.ledger.cumulative_dissipated_j
                - self.load.state.cumulative_dissipated_j(),
            upstream_work_j: drive.upstream_pressure_pa
                * ((aperture.bore_flow_m3_s - drive.body_flow_m3_s) * dt),
            body_work_j: aperture.bore_pressure_pa * (drive.body_flow_m3_s * dt),
            port_work_residual_j: duct.ledger.port_work_step_j
                - aperture.bore_pressure_pa * (aperture.bore_flow_m3_s * dt),
            duct,
        };
        if ![
            frame.stored_energy_j,
            frame.storage_change_j,
            frame.inferred_total_dissipation_j,
            frame.upstream_work_j,
            frame.body_work_j,
            frame.port_work_residual_j,
            frame.balance_residual_j(),
            frame.duct.ledger.cumulative_port_work_j,
            frame.duct.ledger.cumulative_dissipated_j,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(invalid(
                "DG aperture candidate left the finite energy/work domain",
            ));
        }
        self.load
            .propagator
            .sync_nodal(&mut candidate)
            .map_err(AcousticRealizeError::Duct)?;
        if !candidate
            .pressure()
            .iter()
            .chain(candidate.volume_flow())
            .all(|value| value.is_finite())
        {
            return Err(invalid(
                "DG aperture candidate has nonfinite physical nodal fields",
            ));
        }
        self.load.state = candidate;
        self.aperture.accept_frame(trial);
        Ok(frame)
    }
}

#[cfg(test)]
mod tests;
