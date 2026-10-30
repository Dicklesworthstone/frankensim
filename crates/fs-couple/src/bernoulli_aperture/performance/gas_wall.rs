//! Reciprocal finite gas/wall conduction on the original performance clock.
//!
//! Both caloric inventories, the original material card and the complete
//! moving-source/acoustic histories stay in their existing physical owners.
//! Controls prescribe external powers; the shared thermal solver determines
//! interface heat and pays the native gas and wall parameter work.

use super::{AperturePerformance, CoupledAperture, GestureCompileError, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dg::thermal_wall::gas::{
    FiniteGasWallApertureDgStep, FiniteGasWallThermalSolve, admit_initial,
};
use crate::bernoulli_aperture::dg::{ApertureDg, DgApertureFrame};
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, UniformEnthalpyStepInput,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use std::sync::Arc;

/// Owned declaration of the existing material series-conduction interface.
///
/// The retained original wall card supplies solid conductivity; the same
/// thermochemical gas that owns internal energy supplies gas conductivity.
/// Each coefficient is evaluated at its own actual energy midpoint. Effective
/// distances describe unresolved conduction layers, not convection or a
/// resolved temperature profile.
#[derive(Clone, Debug, PartialEq)]
pub struct GasWallConduction {
    /// All original external material conditions, with explicit quantity axes.
    pub context: QueryPoint,
    /// Exact scalar, static thermal-conductivity property and quantity.
    pub conductivity: PropertyKey,
    /// Original source selection or explicit property pin.
    pub selection: MaterialPropertySelection,
    /// Optional existing liquid-mass-fraction coordinate; this route is solid.
    pub liquid_fraction_axis: Option<String>,
    /// Common transfer area in square metres; zero explicitly means isolation.
    pub area_m2: f64,
    /// Positive effective gas conduction distance in metres.
    pub gas_distance_m: f64,
    /// Positive effective solid conduction distance in metres.
    pub body_distance_m: f64,
}

impl GasWallConduction {
    fn layers<'a>(
        &'a self,
        material: &'a MaterialCard,
        temperature_axis: &'a str,
    ) -> MaterialConductionLayers<'a> {
        MaterialConductionLayers {
            material,
            context: &self.context,
            temperature_axis,
            liquid_fraction_axis: self.liquid_fraction_axis.as_deref(),
            conductivity: &self.conductivity,
            selection: self.selection.clone(),
            area_m2: self.area_m2,
            gas_distance_m: self.gas_distance_m,
            body_distance_m: self.body_distance_m,
        }
    }
}

/// Paired signed external powers applied at one physical source sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasWallHeatControl {
    /// Sample in the performance's half-open source window.
    pub sample: u64,
    /// Watts into the closed gas inventory; negative values prescribe cooling.
    pub gas_power_w: f64,
    /// Watts into the finite wall body; negative values prescribe cooling.
    pub wall_power_w: f64,
}

/// External power history and explicit coupled material/energy budgets.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasWallHeatProgram {
    /// Strictly increasing samples beginning at zero. Both powers are held
    /// until the next event. Each pair consumes two shared control/work slots.
    pub controls: Vec<GasWallHeatControl>,
    /// Existing gas-energy, material-work and interface residual admissions.
    pub coupling: FiniteGasWallThermalSolve,
}

pub(super) struct FiniteGasWallPlayback {
    gas: ClosedGasParcel,
    curve: EquilibriumEnthalpyPhaseCurve,
    state: EquilibriumPhaseState,
    material: Arc<MaterialCard>,
    temperature_axis: String,
    conduction: GasWallConduction,
    mass_kg: f64,
    volume_m3: f64,
    surface_area_m2: f64,
    duration_s: f64,
    program: FiniteGasWallHeatProgram,
    gate: Arc<CancelGate>,
    next: usize,
    held_power_w: [f64; 2],
    last: Option<FiniteGasWallApertureDgStep>,
}

impl FiniteGasWallPlayback {
    fn candidate(&self, sample: u64) -> [f64; 2] {
        self.program
            .controls
            .get(self.next)
            .filter(|event| event.sample == sample)
            .map_or(self.held_power_w, |event| {
                [event.gas_power_w, event.wall_power_w]
            })
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

    pub(super) fn assignment_count(&self) -> usize {
        // Multiplication and shared capacity were checked before attachment.
        self.program.controls.len() * 2
    }

    pub(super) fn step(
        &mut self,
        system: &mut ApertureDg,
        drive: TubeDrive,
        force_n: f64,
    ) -> Result<DgApertureFrame, RenderError> {
        let powers = self.candidate(system.aperture().accepted_steps());
        let body = UniformEnthalpyStepInput {
            curve: &self.curve,
            initial: self.state,
            mass_kg: self.mass_kg,
            volume_m3: self.volume_m3,
            surface_area_m2: self.surface_area_m2,
            internal_heat_j: powers[1] * self.duration_s,
            duration_s: self.duration_s,
        };
        let layers = self
            .conduction
            .layers(&self.material, &self.temperature_axis);
        let report = system
            .step_with_finite_gas_wall_heat_controlled(
                &mut self.gas,
                &mut self.state,
                body,
                &layers,
                drive,
                force_n,
                powers[0] * self.duration_s,
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
            self.held_power_w = [event.gas_power_w, event.wall_power_w];
            self.next += 1;
        }
    }
}

impl AperturePerformance {
    /// Evolve both finite inventories through reciprocal material conduction.
    ///
    /// The input declares the original wall enthalpy chart, accepted initial
    /// state and complete fixed main-wall inventory. Its duration must equal
    /// the mechanical clock and its internal heat must be zero; the program
    /// supplies external powers. The original immutable card is retained by Arc
    /// and must own both this chart and every selected wall region.
    ///
    /// One mechanical/acoustic march supplies native solid dissipation. The
    /// joint transaction closes gas energy, body enthalpy, native gas/wall
    /// parameter work and midpoint conduction. Interface heat is credited to
    /// the gas and debited from the body once. Deposited native solid heat is
    /// removed once from exported loss. Material updates retain the original
    /// acoustic coordinates, wall memories, terminal histories and clocks.
    ///
    /// This is one closed, fixed-composition, fixed-volume gas inventory and
    /// one all-solid wall body with fixed geometry and density. It retains the
    /// original first-order frozen acoustic partition. Independent upstream or
    /// branch inventories, computed exterior mean-pressure feedback, expansion,
    /// phase transfer, advection and resolved convection need separate owners.
    /// It cannot be attached alongside another evolving gas/wall program.
    ///
    /// Admission advances no accepted state. A refused physical step preserves
    /// both inventories and all physical histories. Controls advance only after
    /// successful observation; an observation failure retains the existing
    /// renderer's poisoned-callback semantics.
    #[allow(clippy::too_many_arguments)]
    pub fn with_finite_gas_wall_heat(
        mut self,
        gas: ClosedGasParcel,
        body: UniformEnthalpyStepInput<'_>,
        material: Arc<MaterialCard>,
        temperature_axis: String,
        conduction: GasWallConduction,
        program: FiniteGasWallHeatProgram,
        gate: Arc<CancelGate>,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned
            || self.system.aperture().accepted_steps() != 0
            || self.finite_gas_wall.is_some()
            || self.finite_gas.is_some()
            || self.finite_transport.is_some()
            || self.finite_wall.is_some()
        {
            return Err(invalid(
                "joint gas/wall heat must bind once before playback, with one evolving owner per inventory",
            ));
        }
        if body.internal_heat_j != 0.0
            || body.duration_s.to_bits() != self.system.aperture().spec().time_step_s.to_bits()
        {
            return Err(invalid(
                "joint wall input must have the mechanical duration and zero heat before its power program is attached",
            ));
        }
        let (assignments, work) = self.remaining_control_budgets()?;
        let needed = program
            .controls
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("joint gas/wall heat assignment count overflow"))?;
        if program.controls.is_empty()
            || needed > assignments
            || needed as u128 > u128::from(work)
            || program.controls[0].sample != 0
        {
            return Err(invalid(
                "joint gas/wall heat needs an initial pair and shared control/work capacity for both powers",
            ));
        }
        let mut previous = None;
        for event in &program.controls {
            if gate.is_requested() {
                return Err(GestureCompileError::Render(RenderError::Voice(
                    AcousticRealizeError::Cancelled,
                )));
            }
            if event.sample >= self.config.samples
                || previous.is_some_and(|sample| event.sample <= sample)
            {
                return Err(invalid(
                    "joint heat assignments require strictly increasing in-window samples",
                ));
            }
            for power in [event.gas_power_w, event.wall_power_w] {
                let heat = power * body.duration_s;
                if !power.is_finite() || !heat.is_finite() || (power != 0.0 && heat == 0.0) {
                    return Err(invalid(
                        "joint heat assignments require representable signed gas and wall powers",
                    ));
                }
            }
            previous = Some(event.sample);
        }
        let CoupledAperture::Dg(system) = &self.system else {
            return Err(invalid(
                "joint gas/wall heat requires the original source-bound DG aperture owner",
            ));
        };
        let layers = conduction.layers(&material, &temperature_axis);
        admit_initial(system, &gas, body, &layers, program.coupling, &mut || {
            gate.is_requested()
        })
        .map_err(|error| GestureCompileError::Render(RenderError::Voice(error)))?;
        self.finite_gas_wall = Some(FiniteGasWallPlayback {
            gas,
            curve: body.curve.clone(),
            state: body.initial,
            material,
            temperature_axis,
            conduction,
            mass_kg: body.mass_kg,
            volume_m3: body.volume_m3,
            surface_area_m2: body.surface_area_m2,
            duration_s: body.duration_s,
            program,
            gate,
            next: 0,
            held_power_w: [0.0; 2],
            last: None,
        });
        Ok(self)
    }

    /// Actual accepted gas inventory and its original thermochemical source.
    #[must_use]
    pub fn finite_gas_wall_gas(&self) -> Option<&ClosedGasParcel> {
        self.finite_gas_wall.as_ref().map(|joint| &joint.gas)
    }

    /// Accepted finite wall input and original chart, with no next heat imposed.
    #[must_use]
    pub fn finite_gas_wall_body(&self) -> Option<UniformEnthalpyStepInput<'_>> {
        self.finite_gas_wall
            .as_ref()
            .map(FiniteGasWallPlayback::body)
    }

    /// Actual accepted wall enthalpy state.
    #[must_use]
    pub fn finite_gas_wall_state(&self) -> Option<&EquilibriumPhaseState> {
        self.finite_gas_wall.as_ref().map(|joint| &joint.state)
    }

    /// Original immutable wall card used by the mechanical and thermal queries.
    #[must_use]
    pub fn finite_gas_wall_material(&self) -> Option<&MaterialCard> {
        self.finite_gas_wall
            .as_ref()
            .map(|joint| joint.material.as_ref())
    }

    /// Authored interface geometry, source conditions and conductivity key.
    #[must_use]
    pub fn finite_gas_wall_conduction(&self) -> Option<&GasWallConduction> {
        self.finite_gas_wall.as_ref().map(|joint| &joint.conduction)
    }

    /// Original paired external powers and coupling budgets.
    #[must_use]
    pub fn finite_gas_wall_heat_program(&self) -> Option<&FiniteGasWallHeatProgram> {
        self.finite_gas_wall.as_ref().map(|joint| &joint.program)
    }

    /// Shared cancellation gate retained throughout every material iteration.
    #[must_use]
    pub fn finite_gas_wall_cancel_gate(&self) -> Option<&Arc<CancelGate>> {
        self.finite_gas_wall.as_ref().map(|joint| &joint.gate)
    }

    /// Power pairs published with successfully observed physical samples.
    #[must_use]
    pub fn applied_gas_wall_heat_controls(&self) -> &[GasWallHeatControl] {
        self.finite_gas_wall
            .as_ref()
            .map_or(&[], |joint| &joint.program.controls[..joint.next])
    }

    /// Complete coupled energy report for the last accepted physical step.
    #[must_use]
    pub fn last_finite_gas_wall_step(&self) -> Option<&FiniteGasWallApertureDgStep> {
        self.finite_gas_wall
            .as_ref()
            .and_then(|joint| joint.last.as_ref())
    }
}
