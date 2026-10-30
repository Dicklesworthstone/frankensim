//! Scheduled heat in the original finite material wall and moving aperture.
//!
//! The existing aperture/load transaction supplies solved flows, native wall
//! dissipation and signed material work. One source-bound enthalpy chart owns
//! the body temperature. Controls select external power, never a target pitch
//! or temperature, and publish only with their accepted physical sample.

use super::{AperturePerformance, CoupledAperture, GestureCompileError, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dg::thermal_wall::{
    FiniteWallApertureDgStep, FiniteWallThermalSolve, admit_initial,
};
use crate::bernoulli_aperture::dg::{ApertureDg, DgApertureFrame};
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_matdb::MaterialCard;
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, UniformEnthalpyStepInput,
};
use std::sync::Arc;

/// Signed external power supplied to the retained finite wall body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallHeatControl {
    /// Physical source sample in the performance's half-open window.
    pub sample: u64,
    /// Watts into the body. Negative values prescribe external cooling.
    pub power_w: f64,
}

/// External heat controls and explicit material-work/source-heat budgets.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteWallHeatProgram {
    /// Strictly increasing samples starting at zero; the last power is held.
    pub controls: Vec<WallHeatControl>,
    /// Existing scalar energy closure and native material/heat admissions.
    pub coupling: FiniteWallThermalSolve,
}

pub(super) struct FiniteWallPlayback {
    curve: EquilibriumEnthalpyPhaseCurve,
    state: EquilibriumPhaseState,
    material: Arc<MaterialCard>,
    temperature_axis: String,
    mass_kg: f64,
    volume_m3: f64,
    surface_area_m2: f64,
    duration_s: f64,
    pub(super) program: FiniteWallHeatProgram,
    gate: Arc<CancelGate>,
    next: usize,
    held_power_w: f64,
    last: Option<FiniteWallApertureDgStep>,
}

impl FiniteWallPlayback {
    fn candidate(&self, sample: u64) -> f64 {
        self.program
            .controls
            .get(self.next)
            .filter(|event| event.sample == sample)
            .map_or(self.held_power_w, |event| event.power_w)
    }

    fn body(&self) -> UniformEnthalpyStepInput<'_> {
        UniformEnthalpyStepInput {
            curve: &self.curve,
            initial: self.state,
            mass_kg: self.mass_kg,
            volume_m3: self.volume_m3,
            surface_area_m2: self.surface_area_m2,
            internal_heat_j: 0.0,
            duration_s: self.duration_s,
        }
    }

    pub(super) fn step(
        &mut self,
        system: &mut ApertureDg,
        drive: TubeDrive,
        force_n: f64,
    ) -> Result<DgApertureFrame, RenderError> {
        let heat_j = self.candidate(system.aperture().accepted_steps()) * self.duration_s;
        let body = UniformEnthalpyStepInput {
            curve: &self.curve,
            initial: self.state,
            mass_kg: self.mass_kg,
            volume_m3: self.volume_m3,
            surface_area_m2: self.surface_area_m2,
            internal_heat_j: heat_j,
            duration_s: self.duration_s,
        };
        let report = system
            .step_with_finite_wall_heat_controlled(
                &mut self.state,
                body,
                &self.material,
                &self.temperature_axis,
                drive,
                force_n,
                self.program.coupling,
                &mut || self.gate.is_requested(),
            )
            .map_err(RenderError::Voice)?;
        let frame = report.coupled.clone();
        self.last = Some(report);
        Ok(frame)
    }

    pub(super) fn accept(&mut self, sample: u64) {
        if let Some(event) = self
            .program
            .controls
            .get(self.next)
            .filter(|event| event.sample == sample)
        {
            self.held_power_w = event.power_w;
            self.next += 1;
        }
    }
}

impl AperturePerformance {
    /// Schedule physical heat into the original source-bound DG wall body.
    ///
    /// The supplied input declares the accepted initial enthalpy, fixed
    /// inventory and original chart. Its duration must equal the mechanical
    /// clock and its internal_heat_j must be zero: the program alone supplies
    /// external power during playback. The chart is retained by value; the
    /// immutable original material card is retained through its Arc.
    ///
    /// Uniform and regional wall bindings keep their original source queries,
    /// claims, thicknesses, restraints and native histories. Regions share this
    /// one temperature and caloric chart. The body must contain their complete
    /// modeled material inventory and remain all-solid at constant density.
    /// Each step closes body enthalpy against native signed wall work. Optional
    /// Kelvin/Maxwell heat comes from the original solved-flow source integral;
    /// gas, numerical, mechanical and terminal losses are not inferred as heat.
    ///
    /// The body owns the main bore wall; finite side branches retain their
    /// independent materials. The complete downstream terminal/branch load and
    /// its exterior observation remain active. This fixed-gas route does not add gas/body
    /// conduction, independent upstream heating, expansion or phase transfer.
    /// A finite-gas heat program requires a joint thermal transaction and cannot
    /// be attached simultaneously.
    ///
    /// Admission never advances physics. A refused physical sample leaves its
    /// body, aperture and complete load unchanged; the existing renderer then
    /// poisons the callback. Earlier completed samples remain inspectable.
    pub fn with_finite_wall_heat(
        mut self,
        body: UniformEnthalpyStepInput<'_>,
        material: Arc<MaterialCard>,
        temperature_axis: String,
        program: FiniteWallHeatProgram,
        gate: Arc<CancelGate>,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned
            || self.system.aperture().accepted_steps() != 0
            || self.finite_wall.is_some()
            || self.finite_gas.is_some()
        {
            return Err(invalid(
                "finite wall heat must bind once before playback; combined finite gas/body heat needs a joint transaction",
            ));
        }
        if body.internal_heat_j != 0.0
            || body.duration_s.to_bits() != self.system.aperture().spec().time_step_s.to_bits()
        {
            return Err(invalid(
                "finite wall input must have the mechanical duration and zero heat before its power program is attached",
            ));
        }
        let (assignments, work) = self.remaining_control_budgets()?;
        if program.controls.is_empty()
            || program.controls.len() > assignments
            || program.controls.len() as u128 > u128::from(work)
            || program.controls[0].sample != 0
        {
            return Err(invalid(
                "finite wall heat needs an initial assignment and shared control/work capacity",
            ));
        }
        let mut previous = None;
        for event in &program.controls {
            if gate.is_requested() {
                return Err(GestureCompileError::Render(RenderError::Voice(
                    AcousticRealizeError::Cancelled,
                )));
            }
            let heat_j = event.power_w * body.duration_s;
            if event.sample >= self.config.samples
                || previous.is_some_and(|sample| event.sample <= sample)
                || !event.power_w.is_finite()
                || !heat_j.is_finite()
                || (event.power_w != 0.0 && heat_j == 0.0)
            {
                return Err(invalid(
                    "wall heat assignments require increasing in-window samples and representable signed power",
                ));
            }
            previous = Some(event.sample);
        }
        let CoupledAperture::Dg(system) = &self.system else {
            return Err(invalid(
                "finite wall heating requires the original source-bound DG aperture owner",
            ));
        };
        admit_initial(
            system,
            body,
            &material,
            &temperature_axis,
            program.coupling,
            &mut || gate.is_requested(),
        )
        .map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
        self.finite_wall = Some(FiniteWallPlayback {
            curve: body.curve.clone(),
            state: body.initial,
            material,
            temperature_axis,
            mass_kg: body.mass_kg,
            volume_m3: body.volume_m3,
            surface_area_m2: body.surface_area_m2,
            duration_s: body.duration_s,
            program,
            gate,
            next: 0,
            held_power_w: 0.0,
            last: None,
        });
        Ok(self)
    }

    /// Accepted finite body and its original chart, with no prescribed next heat.
    #[must_use]
    pub fn finite_wall_body(&self) -> Option<UniformEnthalpyStepInput<'_>> {
        self.finite_wall.as_ref().map(FiniteWallPlayback::body)
    }

    /// Actual accepted body state, shared by every selected regional wall.
    #[must_use]
    pub fn finite_wall_state(&self) -> Option<&EquilibriumPhaseState> {
        self.finite_wall.as_ref().map(|wall| &wall.state)
    }

    /// Original immutable source card used for every temperature refresh.
    #[must_use]
    pub fn finite_wall_material(&self) -> Option<&MaterialCard> {
        self.finite_wall.as_ref().map(|wall| wall.material.as_ref())
    }

    /// Original authored heat powers and coupling admissions.
    #[must_use]
    pub fn finite_wall_heat_program(&self) -> Option<&FiniteWallHeatProgram> {
        self.finite_wall.as_ref().map(|wall| &wall.program)
    }

    /// Cancellation gate retained throughout the physical feedback solve.
    #[must_use]
    pub fn finite_wall_cancel_gate(&self) -> Option<&Arc<CancelGate>> {
        self.finite_wall.as_ref().map(|wall| &wall.gate)
    }

    /// Heat controls published with successfully observed physical samples.
    #[must_use]
    pub fn applied_wall_heat_controls(&self) -> &[WallHeatControl] {
        self.finite_wall
            .as_ref()
            .map_or(&[], |wall| &wall.program.controls[..wall.next])
    }

    /// Native mechanical, acoustic, wall and finite-body account of the last step.
    #[must_use]
    pub fn last_finite_wall_step(&self) -> Option<&FiniteWallApertureDgStep> {
        self.finite_wall
            .as_ref()
            .and_then(|wall| wall.last.as_ref())
    }
}
