//! Several structural modes coupled through ONE slit and pressure/flow port.
//!
//! The gap is H + g.q, pressure forces are -A*dp and swept flow is -A.v.
//! A and g are separate signed projections; neither is a fitted sound gain.
//! Contact and continuous Bernoulli flow remain in their original owners. The
//! scalar solve generalizes the existing independent-blade response condensation.
//! This is a fixed linear modal basis and a uniform effective slit. An optional
//! full passive viscosity matrix retains cross-mode loss in the SAME solve.
//! ExactZoh advances the original structural pencil with held pressure/contact
//! forces and integrated material dissipation; nonlinear slit/contact time error
//! still needs refinement. Proportional source-material Maxwell spectra retain
//! every mode and their actual histories in that same exact solve.

use super::dynamic::{ApertureDrive, ApertureFrame, ApertureState};
use super::paired::{BladePropagator, BladeState, PairedApertureLimits};
use super::uniform::{self, MAX_MODES};
use crate::acoustic_realize::AcousticRealizeError;
use crate::unilateral_contact::distributed::ApertureContactStep;
use fs_dcontact::{ContactStorage, Obstacle};
use fs_phs::Storage;
use fs_phs::{
    NewtonSite, SolveBudget, SolveCertificate, admit_newton_iterate, arithmetic_rounding_floor,
};
mod exact;
/// Source-bound nodal clearance, contact and physical slit flow.
pub mod closure;
mod spatial;
/// Original regional material plates and same-basis spatial ports.
pub mod material;
/// Prescribed spatial pressure loads on the original material modal basis.
pub mod pressure;
/// Original wetted face geometry in the coupled acoustic force/volume port.
pub mod wetted;
/// Source-bound scalar, regional and full-tensor Maxwell material memory.
pub mod relaxation;
mod viscosity;
/// Complete tetrahedral material images in the original coupled modal owner.
pub mod solid;

/// One fixed-basis mode, in a physical metre-valued coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ApertureMode {
    /// Positive generalized mass [kg].
    pub mass_kg: f64,
    /// Positive generalized restoring stiffness [N/m]; the equilibrium
    /// stiffness when source-bound Maxwell branches are selected.
    pub stiffness_n_m: f64,
    /// Explicit nonnegative diagonal viscous loss [N s/m].
    pub damping_n_s_m: f64,
    /// Signed pressure virtual-work area and swept-volume projection [m²].
    pub pressure_area_m2: f64,
    /// Signed mean-slit displacement per coordinate metre [dimensionless].
    pub gap_projection: f64,
    /// Two independent physical force ports: generalized force is this row
    /// dotted with the supplied newtons. An all-zero column is an absent port.
    pub force_projection: [f64; 2],
}

/// Immutable physical basis, fluid load, clock and existing explicit limits.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalApertureSpec {
    /// Positive effective rest gap [m].
    pub rest_opening_m: f64,
    /// Single shared slit width [m], not multiplied by the mode count.
    pub width_m: f64,
    /// 1..=32 independent modes. Missing modes are NOT inferred or rescaled.
    pub modes: Vec<ApertureMode>,
    /// Actual inlet gas density [kg/m³].
    pub density_kg_m3: f64,
    /// Characteristic pressure/flow load [Pa s/m³].
    pub impedance_pa_s_m3: f64,
    /// Fixed mechanical step [s].
    pub time_step_s: f64,
    /// Total accepted samples, including continuation; at most 2^53.
    pub max_steps: u64,
    /// Bounds apply to every mode; velocity tolerance applies to the total gap.
    pub limits: PairedApertureLimits,
}

/// Completely accepted mechanics and original pressure/flow work account.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalApertureFrame {
    /// Common slit, fluid and complete mechanical/contact/material energy observation.
    pub aperture: ApertureFrame,
    /// Independent endpoint coordinates in the original basis order.
    pub modes: Vec<BladeState>,
    /// Same-sample pressure difference across the source [Pa].
    pub pressure_drop_pa: f64,
    /// Sum of physical normal contact forces [N]. A uniform slit uses g_i*F;
    /// spatial closure pulls each point reaction into its own original mode row.
    pub contact_force_n: f64,
    /// Independent force-balance residuals, in the original coordinates [N].
    pub momentum_residual_n: Vec<f64>,
    /// Observed work closure and its declared energy budget [J].
    pub work_certificate: SolveCertificate,
    /// Observed flow closure and its declared volume-flow budget [m^3/s].
    pub flow_certificate: SolveCertificate,
    /// Additional spatial pressure ports, separate from the uniform bore load.
    /// Their work is already included in aperture.mechanical_work_j.
    pub spatial_pressure: [pressure::MaterialPressureFrame; 2],
    port_transfers: [crate::port_adoption::PortTransfer; 6],
    memory: Option<Vec<Vec<f64>>>,
}

/// A complete modal state; no acoustic network or filter is duplicated here.
pub struct ModalAperture {
    spec: ModalApertureSpec,
    state: Vec<BladeState>,
    obstacle: Obstacle,
    contact_storage: ContactStorage,
    closure: Option<closure::ModalPlateClosure>,
    solid_closure: Option<solid::SolidSurfaceClosureSpec>,
    samples: u64,
    material_plates: Vec<material::MaterialModalPlate>,
    material_solid: Option<solid::MaterialModalSolid>,
    hereditary_solid: Option<solid::HereditaryModalSolid>,
    pressure_ports: [Option<pressure::MaterialPressurePort>; 2],
    material_gas: Option<fs_material::gas::GasState>,
    viscous_slit: Option<super::viscous::ViscousSlit>,
    viscosity: Option<viscosity::ApertureViscosity>,
    exact: Option<exact::ExactModes>,
    relaxation: Option<relaxation::ModalRelaxation>,
    ports: crate::port_adoption::ModulePorts,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn nonlinear(e: impl core::fmt::Display) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(e.to_string())
}
struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 {
        0.0
    }
    fn gradient(&self, _: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }
}

impl ModalAperture {
    /// Admit a supplied physical modal basis and ONE effective gap contact.
    /// There is no named-instrument preset or inferred damping. Modal sign and
    /// scaling belong to the supplied projections, including the force ports.
    /// A successful construction does not validate omitted modes or mean-gap
    /// closure. Select full viscosity explicitly with `with_viscosity`;
    /// source-bound proportional hereditary laws use `with_material_relaxation`
    /// after choosing the exact structural propagator.
    pub fn new(
        spec: ModalApertureSpec,
        state: Vec<BladeState>,
        obstacle: Obstacle,
    ) -> Result<Self, AcousticRealizeError> {
        Self::new_with_propagator(spec, state, obstacle, BladePropagator::Midpoint)
    }

    /// Select the original midpoint image or exact zero-order-held structural
    /// propagation, using the same pressure, swept-volume and contact solve.
    /// ExactZoh retains the physical M/K/C and free-mode decay without frequency
    /// warping. Pressure, contact and supplied forces are held over each sample;
    /// nonlinear/contact accuracy and sample aliasing still require clock refinement.
    /// Full viscosity selected subsequently is included in the SAME exact block.
    pub fn new_with_propagator(
        spec: ModalApertureSpec,
        state: Vec<BladeState>,
        obstacle: Obstacle,
        propagator: BladePropagator,
    ) -> Result<Self, AcousticRealizeError> {
        let l = spec.limits;
        if !(1..=MAX_MODES).contains(&spec.modes.len())
            || state.len() != spec.modes.len()
            || ![
                spec.rest_opening_m,
                spec.width_m,
                spec.density_kg_m3,
                spec.impedance_pa_s_m3,
                spec.time_step_s,
                l.maximum_abs_displacement_m,
                l.maximum_abs_velocity_m_s,
                l.maximum_energy_j,
                l.velocity_tolerance_m_s,
                l.relative_tolerance,
            ]
            .iter()
            .all(|x| x.is_finite() && *x > 0.0)
            || l.relative_tolerance >= 1.0 / 128.0
            || !(1..=512).contains(&l.max_evaluations)
            || spec.max_steps == 0
            || spec.max_steps > (1_u64 << 53)
            || !(spec.time_step_s * spec.max_steps as f64).is_finite()
            || obstacle.provenance().trim().is_empty()
        {
            return Err(invalid(
                "modal aperture requires bounded modes, matched states and finite physical limits",
            ));
        }
        for mode in &spec.modes {
            if ![mode.mass_kg, mode.stiffness_n_m]
                .iter()
                .all(|x| x.is_finite() && *x > 0.0)
                || !mode.damping_n_s_m.is_finite()
                || mode.damping_n_s_m < 0.0
                || [
                    mode.pressure_area_m2,
                    mode.gap_projection,
                    mode.force_projection[0],
                    mode.force_projection[1],
                ]
                .iter()
                .any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "modal aperture has invalid mechanics or signed physical projections",
                ));
            }
        }
        if spec.modes.iter().all(|m| m.gap_projection == 0.0) {
            return Err(invalid(
                "modal aperture basis does not move its declared slit",
            ));
        }
        let gap = spec.rest_opening_m
            + spec
                .modes
                .iter()
                .zip(&state)
                .map(|(m, s)| m.gap_projection * s.displacement_m)
                .sum::<f64>();
        ApertureContactStep::new(&obstacle, gap).map_err(nonlinear)?;
        let contact_storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![obstacle.clone()])
            .map_err(nonlinear)?;
        let ports = crate::port_adoption::ModulePorts::new(
            "modal-aperture",
            vec![
                ("upstream", fs_ports::PortKind::FluidPressureFlux),
                ("downstream", fs_ports::PortKind::FluidPressureFlux),
                ("plate/0", fs_ports::PortKind::MechanicalForceVelocity),
                ("plate/1", fs_ports::PortKind::MechanicalForceVelocity),
            ],
        );
        let exact = match propagator {
            BladePropagator::Midpoint => None,
            BladePropagator::ExactZoh => Some(exact::ExactModes::new(&spec, None)?),
        };
        let model = Self {
            spec,
            state,
            obstacle,
            contact_storage,
            closure: None,
            solid_closure: None,
            samples: 0,
            material_plates: Vec::new(),
            material_solid: None,
            hereditary_solid: None,
            pressure_ports: [None, None],
            material_gas: None,
            viscous_slit: None,
            viscosity: None,
            exact,
            relaxation: None,
            ports,
        };
        model.admit_linear_clock()?;
        model.validate_state(&model.state)?;
        Ok(model)
    }
    /// Original physical inputs and explicit approximation limits.
    #[must_use]
    pub const fn spec(&self) -> &ModalApertureSpec {
        &self.spec
    }
    /// Actual independent physical coordinates, not a single equivalent reed.
    #[must_use]
    pub fn modal_states(&self) -> &[BladeState] {
        &self.state
    }
    /// Accepted mechanical samples only.
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 {
        self.samples
    }
    /// Selected linear structural propagator; no material or modal state changes.
    #[must_use]
    pub const fn blade_propagator(&self) -> BladePropagator {
        if self.exact.is_some() {
            BladePropagator::ExactZoh
        } else {
            BladePropagator::Midpoint
        }
    }

    /// Select the common channel law before stepping. Every mode sees the SAME
    /// pressure drop and its own conjugate swept area; adding modes never adds
    /// parallel slits or multiplies the physical flow width.
    /// A source-bound law retains its full gas for subsequent network admission;
    /// it must agree with any gas already supplied by the material construction.
    pub fn set_viscous_slit(
        &mut self,
        law: super::viscous::ViscousSlit,
    ) -> Result<(), AcousticRealizeError> {
        if self.samples != 0 {
            return Err(invalid(
                "modal viscous slit must be selected before the first accepted sample",
            ));
        }
        law.evaluate(
            self.spec.width_m,
            self.gap(&self.state).max(0.0),
            self.spec.density_kg_m3,
            0.0,
        )?;
        // Preserve the complete source binding for connected network admission.
        // Density alone cannot identify viscosity, heat capacity or temperature.
        if let Some(gas) = law.gas().copied() {
            if self.material_gas.is_some_and(|retained| retained != gas) {
                return Err(invalid(
                    "modal viscous slit must use the complete retained source gas",
                ));
            }
            self.material_gas = Some(gas);
        }
        self.viscous_slit = Some(law);
        Ok(())
    }

    /// Derive channel viscosity from the retained material gas and supplied
    /// streamwise geometry. Missing source gas is refused, never replaced by air.
    pub fn set_material_viscous_slit(&mut self, length_m: f64) -> Result<(), AcousticRealizeError> {
        let gas = self
            .material_gas
            .as_ref()
            .ok_or_else(|| invalid("modal channel requires its retained material gas"))?;
        self.set_viscous_slit(super::viscous::ViscousSlit::from_gas(gas, length_m)?)
    }

    /// Original explicit or source-bound channel law.
    #[must_use]
    pub const fn viscous_slit(&self) -> Option<super::viscous::ViscousSlit> {
        self.viscous_slit
    }

    pub(super) fn rebind_exact_gas_load(
        &mut self, gas: &fs_material::gas::GasState, impedance_pa_s_m3: f64,
    ) -> Result<(), AcousticRealizeError> {
        gas.validate_coherence().map_err(|_| invalid("live modal aperture requires one coherent source gas"))?;
        if self.exact.is_none() || !impedance_pa_s_m3.is_finite() || impedance_pa_s_m3 <= 0.0 {
            return Err(invalid("live modal gas refresh requires ExactZoh mechanics and positive finite impedance"));
        }
        let law = self.viscous_slit.map(|law| law.rebind_source_gas(gas)).transpose()?;
        if let Some(law) = law {
            match &self.closure {
                Some(closure) => {
                    let mut coordinates = [0.0; MAX_MODES];
                    for (q, state) in coordinates.iter_mut().zip(&self.state) {
                        *q = state.displacement_m;
                    }
                    closure.viscous_flow(law, &coordinates[..self.state.len()], gas.density, 0.0)?;
                }
                None => { law.evaluate(self.spec.width_m, self.gap(&self.state).max(0.0), gas.density, 0.0)?; }
            }
        }
        // A gas intervention changes no structural mode, material relaxation,
        // contact coordinate or observation port and advances no sample clock.
        self.spec.density_kg_m3 = gas.density;
        self.spec.impedance_pa_s_m3 = impedance_pa_s_m3;
        self.viscous_slit = law;
        self.material_gas = Some(*gas);
        Ok(())
    }
    fn gap(&self, state: &[BladeState]) -> f64 {
        if let Some(closure) = &self.closure {
            return closure.state_gap(state);
        }
        self.spec.rest_opening_m
            + self
                .spec
                .modes
                .iter()
                .zip(state)
                .map(|(m, s)| m.gap_projection * s.displacement_m)
                .sum::<f64>()
    }
    fn energy(&self, state: &[BladeState]) -> f64 {
        self.closure.as_ref().map_or_else(
            || self.contact_storage.hamiltonian(&[self.gap(state), 0.0]),
            |closure| closure.energy(state))
            + self
                .spec
                .modes
                .iter()
                .zip(state)
                .map(|(m, s)| {
                    0.5 * m.mass_kg * s.velocity_m_s * s.velocity_m_s
                        + 0.5 * m.stiffness_n_m * s.displacement_m * s.displacement_m
                })
                .sum::<f64>()
    }
    /// Mechanical, contact and selected material-memory energy [J], evaluated
    /// from the accepted state rather than accumulated input work.
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.energy(&self.state)
            + self
                .relaxation
                .as_ref()
                .map_or(0.0, |r| r.energy(&self.state, &r.states))
    }
    fn validate_state(&self, state: &[BladeState]) -> Result<(), AcousticRealizeError> {
        let l = self.spec.limits;
        if state.iter().any(|s| {
            !s.displacement_m.is_finite()
                || !s.velocity_m_s.is_finite()
                || s.displacement_m.abs() > l.maximum_abs_displacement_m
                || s.velocity_m_s.abs() > l.maximum_abs_velocity_m_s
        }) {
            return Err(invalid(
                "modal aperture state exceeds its physical coordinate bounds",
            ));
        }
        for plate in &self.material_plates {
            plate.validate_state(state, self.closure.is_none())?;
        }
        if let Some(solid) = &self.material_solid {
            solid.validate_state(state)?;
        }
        if let Some(solid) = &self.hereditary_solid {
            solid.validate_state(state)?;
        }
        if let Some(closure) = &self.closure {
            closure.validate_state(state)?;
        }
        let energy = self.energy(state);
        if !self.gap(state).is_finite()
            || !energy.is_finite()
            || energy < 0.0
            || energy > l.maximum_energy_j
        {
            return Err(invalid(
                "modal aperture gap or stored energy exceeds its physical domain",
            ));
        }
        Ok(())
    }
    /// Check finite forces and reject an input addressed to an absent force port.
    pub fn validate_forces(&self, forces: [f64; 2]) -> Result<(), AcousticRealizeError> {
        for i in 0..2 {
            if !forces[i].is_finite()
                || (forces[i] != 0.0
                    && self.spec.modes.iter().all(|m| m.force_projection[i] == 0.0))
            {
                return Err(invalid(
                    "modal force requires finite newtons and an existing physical projection",
                ));
            }
        }
        Ok(())
    }
    /// Extend a finite horizon without changing any physical coordinate.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        if total <= self.spec.max_steps
            || total > (1_u64 << 53)
            || !(total as f64 * self.spec.time_step_s).is_finite()
        {
            return Err(invalid(
                "modal aperture horizon extension must increase a finite sample budget",
            ));
        }
        self.spec.max_steps = total;
        Ok(())
    }
    /// Advance pressure and all modes through their shared contact and flow solve.
    pub fn step(
        &mut self,
        drive: ApertureDrive,
        forces: [f64; 2],
    ) -> Result<ModalApertureFrame, AcousticRealizeError> {
        let frame = self.preview(drive, forces)?;
        self.accept(&frame);
        Ok(frame)
    }
    pub(super) fn accept(&mut self, frame: &ModalApertureFrame) {
        self.ports.accept(
            crate::port_adoption::StepAccount {
                tick: frame.aperture.step,
                duration_s: self.spec.time_step_s,
                storage_change_j: Some(frame.aperture.storage_change_j),
                loss_j: Some(frame.aperture.dissipated_energy_j),
                included_work_j: 0.0,
                entropy_change_j_k: 0.0,
                tolerance_j: Some(frame.work_certificate.effective_budget),
                unavailable: None,
            },
            frame.port_transfers.into_iter().take(self.ports.schemas.len()),
        );
        self.state.copy_from_slice(&frame.modes);
        if let (Some(relaxation), Some(states)) = (&mut self.relaxation, &frame.memory) {
            relaxation.states.clone_from(states);
        }
        self.samples = frame.aperture.step;
    }
    pub(super) fn preview(
        &self,
        drive: ApertureDrive,
        forces: [f64; 2],
    ) -> Result<ModalApertureFrame, AcousticRealizeError> {
        self.preview_with_pressures(drive, forces, [0.0; 2])
    }
    pub(super) fn preview_with_pressures(
        &self,
        drive: ApertureDrive,
        forces: [f64; 2],
        pressures_pa: [f64; 2],
    ) -> Result<ModalApertureFrame, AcousticRealizeError> {
        self.admit_linear_clock()?;
        self.validate_forces(forces)?;
        self.validate_pressures(pressures_pa)?;
        if self.samples >= self.spec.max_steps
            || ![
                drive.upstream_pressure_pa,
                drive.incoming_pressure_pa,
                drive.body_flow_m3_s,
            ]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(invalid(
                "modal aperture needs finite supplies and an unexhausted horizon",
            ));
        }
        let n = self.state.len();
        let dt = self.spec.time_step_s;
        let h = 0.5 * dt;
        let (mut free, mut mobility, mut area, mut shape, mut external) = (
            [0.0; MAX_MODES],
            [0.0; MAX_MODES],
            [0.0; MAX_MODES],
            [0.0; MAX_MODES],
            [0.0; MAX_MODES],
        );
        for (i, (m, s)) in self.spec.modes.iter().zip(&self.state).enumerate() {
            external[i] = m.force_projection[0] * forces[0] + m.force_projection[1] * forces[1];
            for (port, pressure) in self.pressure_ports.iter().zip(pressures_pa) {
                if let Some(port) = port {
                    // Keep the original arithmetic when this additional source
                    // is absent or released; the uniform pressure face is intact.
                    if pressure != 0.0 {
                        external[i] += port.modal_area_m2()[i] * pressure;
                    }
                }
            }
            if !external[i].is_finite() {
                return Err(invalid("material pressure traction overflowed its modal force"));
            }
            area[i] = m.pressure_area_m2;
            shape[i] = m.gap_projection;
            if self.exact.is_some() {
                continue;
            }
            let d = m.mass_kg + h * m.damping_n_s_m + h * h * m.stiffness_n_m;
            mobility[i] = h / d;
            let rhs = m.mass_kg * s.velocity_m_s - h * m.stiffness_n_m * s.displacement_m
                + h * external[i];
            free[i] = if self.viscosity.as_ref().is_some_and(|c| c.is_active()) {
                rhs
            } else {
                rhs / d
            };
            if !d.is_finite()
                || d <= 0.0
                || !mobility[i].is_finite()
                || mobility[i] <= 0.0
                || !free[i].is_finite()
            {
                return Err(invalid("modal midpoint response is not representable"));
            }
        }
        let gap = self.gap(&self.state);
        let z = self.spec.impedance_pa_s_m3;
        let trial = if let Some(closure) = &self.closure {
            if let Some(exact) = &self.exact {
                exact.free_velocity(
                    &self.state, &external[..n], dt,
                    self.relaxation.as_ref().map(|r| r.states.as_slice()), &mut free[..n],
                )?;
            } else if let Some(c) = self.viscosity.as_ref().filter(|c| c.is_active()) {
                c.solve_free(&mut free[..n])?;
            }
            self.solve_profiled(drive, &free[..n], &mobility[..n], &area[..n], closure)?
        } else if let Some(exact) = &self.exact {
            exact.free_velocity(
                &self.state,
                &external[..n],
                dt,
                self.relaxation.as_ref().map(|r| r.states.as_slice()),
                &mut free[..n],
            )?;
            uniform::solve_response(
                &free[..n],
                exact.pressure_response(),
                exact.contact_response(),
                &area[..n],
                &shape[..n],
                gap,
                dt,
                self.spec.width_m,
                self.spec.density_kg_m3,
                z,
                drive,
                self.spec.limits,
                &self.obstacle,
                self.viscous_slit,
            )?
        } else if let Some(c) = self.viscosity.as_ref().filter(|c| c.is_active()) {
            c.solve_free(&mut free[..n])?;
            uniform::solve_response(
                &free[..n],
                c.pressure_response(),
                c.contact_response(),
                &area[..n],
                &shape[..n],
                gap,
                dt,
                self.spec.width_m,
                self.spec.density_kg_m3,
                z,
                drive,
                self.spec.limits,
                &self.obstacle,
                self.viscous_slit,
            )?
        } else {
            uniform::solve(
                &free[..n],
                &mobility[..n],
                &area[..n],
                &shape[..n],
                gap,
                dt,
                self.spec.width_m,
                self.spec.density_kg_m3,
                z,
                drive,
                self.spec.limits,
                &self.obstacle,
                self.viscous_slit,
            )?
        };
        let mut next: Vec<_> = self
            .state
            .iter()
            .zip(&trial.velocity)
            .map(|(s, v)| BladeState {
                displacement_m: s.displacement_m + dt * v,
                velocity_m_s: 2.0 * v - s.velocity_m_s,
            })
            .collect();
        let next_gap = self.gap(&next);
        let vm: f64 = shape[..n]
            .iter()
            .zip(&trial.velocity)
            .map(|(g, v)| g * v)
            .sum();
        let contact = if let Some(closure) = &self.closure {
            let before = closure::displacements(&self.state);
            let after = closure::displacements(&next);
            closure.response(&before[..n], &after[..n], &trial.velocity[..n])?
        } else {
            let contact = ApertureContactStep::new(&self.obstacle, gap)
                .and_then(|c| c.response(next_gap, vm)).map_err(nonlinear)?;
            let mut forces = [0.0; MAX_MODES];
            for i in 0..n { forces[i] = shape[i] * contact.force; }
            closure::ClosureResponse { forces_n: forces, normal_force_n: contact.force,
                dissipated_power_w: contact.dissipated_power }
        };
        let endpoint = self.exact.as_ref().map(|exact| {
            let held: Vec<_> = (0..n)
                .map(|i| external[i] - area[i] * trial.dp + contact.forces_n[i])
                .collect();
            exact.endpoint(
                &self.state,
                &held,
                &trial.velocity[..n],
                dt,
                self.relaxation.as_ref().map(|r| r.states.as_slice()),
            )
        });
        if let Some(exact) = &endpoint {
            next.copy_from_slice(&exact.state);
        }
        self.validate_state(&next)?;
        if let (Some(relaxation), Some(endpoint)) = (&self.relaxation, &endpoint) {
            relaxation.validate_material_state(&next, &endpoint.memory)?;
        }
        let mut middle = [0.0; MAX_MODES];
        for i in 0..n {
            middle[i] = f64::midpoint(self.state[i].displacement_m, next[i].displacement_m);
        }
        let opening = match &self.closure {
            Some(closure) => closure.open_area_m2(&middle[..n])? / self.spec.width_m,
            None => f64::midpoint(gap, next_gap).max(0.0),
        };
        let jet = match (self.viscous_slit, &self.closure) {
            (Some(law), Some(closure)) => closure.viscous_flow(
                law, &middle[..n], self.spec.density_kg_m3, trial.dp)?.flow_m3_s,
            (Some(law), None) => {
                law.evaluate(
                    self.spec.width_m,
                    opening,
                    self.spec.density_kg_m3,
                    trial.dp,
                )?
                .flow_m3_s
            }
            (None, _) => super::moving::volume_flow(
                self.spec.width_m,
                opening,
                trial.dp,
                self.spec.density_kg_m3,
            ),
        };
        let swept = -area[..n]
            .iter()
            .zip(&trial.velocity)
            .map(|(a, v)| a * v)
            .sum::<f64>();
        let no_drop = drive.upstream_pressure_pa - drive.incoming_pressure_pa;
        let load_pressure = z * (jet + swept + drive.body_flow_m3_s);
        let outgoing = if no_drop.abs().max(trial.dp.abs())
            <= drive.incoming_pressure_pa.abs().max(load_pressure.abs())
        {
            no_drop - trial.dp
        } else {
            drive.incoming_pressure_pa + load_pressure
        };
        let wave = (outgoing - drive.incoming_pressure_pa) / z;
        let mut viscous_force = [0.0; MAX_MODES];
        let viscous = if let Some(c) = self.viscosity.as_ref().filter(|c| c.is_active()) {
            c.force_and_power(&trial.velocity[..n], &mut viscous_force[..n])?
        } else {
            self.spec
                .modes
                .iter()
                .zip(&trial.velocity)
                .map(|(m, v)| m.damping_n_s_m * v * v)
                .sum()
        };
        let mechanical: f64 = external[..n]
            .iter()
            .zip(&trial.velocity)
            .map(|(f, v)| dt * f * v)
            .sum();
        let relaxation_energy = self.relaxation.as_ref().map_or(0.0, |r| {
            r.energy(
                &next,
                &endpoint
                    .as_ref()
                    .expect("material memory has an exact owner")
                    .memory,
            )
        });
        let energy = self.energy(&next) + relaxation_energy;
        let relaxation_loss = endpoint.as_ref().map_or(0.0, |e| e.relaxation_loss_j);
        if !relaxation_energy.is_finite()
            || !relaxation_loss.is_finite()
            || !energy.is_finite()
            || energy > self.spec.limits.maximum_energy_j
        {
            return Err(invalid(
                "modal material memory exceeds the finite energy domain",
            ));
        }
        let aperture = ApertureFrame {
            step: self.samples + 1,
            time_s: (self.samples + 1) as f64 * dt,
            state: ApertureState {
                opening_m: next_gap,
                opening_velocity_m_s: match &self.closure {
                    Some(closure) => closure.state_gap_velocity(&next)?,
                    None => shape[..n].iter().zip(&next)
                        .map(|(g, s)| g * s.velocity_m_s).sum(),
                },
            },
            midpoint_opening_m: opening,
            outgoing_pressure_pa: outgoing,
            bore_pressure_pa: outgoing + drive.incoming_pressure_pa,
            jet_flow_m3_s: jet,
            swept_flow_m3_s: swept,
            bore_flow_m3_s: wave,
            flow_residual_m3_s: jet + swept + drive.body_flow_m3_s - wave,
            stored_energy_j: energy,
            storage_change_j: energy - self.stored_energy_j(),
            dissipated_energy_j: match &endpoint {
                Some(exact) => {
                    dt * (trial.dp * jet + contact.dissipated_power_w)
                        + exact.viscous_loss_j
                        + relaxation_loss
                }
                None => dt * (trial.dp * jet + viscous + contact.dissipated_power_w),
            },
            pressure_work_j: dt * trial.dp * (wave - drive.body_flow_m3_s),
            mechanical_work_j: mechanical,
            relaxation_energy_j: relaxation_energy,
            relaxation_loss_j: relaxation_loss,
        };
        let momentum: Vec<_> = self
            .spec
            .modes
            .iter()
            .enumerate()
            .map(|(i, m)| {
                let mean_displacement = endpoint.as_ref().map_or_else(
                    || f64::midpoint(self.state[i].displacement_m, next[i].displacement_m),
                    |exact| exact.mean_displacement[i],
                );
                let residual = m.mass_kg * (next[i].velocity_m_s - self.state[i].velocity_m_s) / dt
                    + m.stiffness_n_m * mean_displacement
                    + m.damping_n_s_m * trial.velocity[i]
                    + endpoint.as_ref().map_or(0.0, |e| e.relaxation_force_n[i])
                    + area[i] * trial.dp
                    - contact.forces_n[i]
                    - external[i];
                if self.viscosity.as_ref().is_some_and(|c| c.is_active()) {
                    residual + viscous_force[i]
                } else {
                    residual
                }
            })
            .collect();
        let flow_scale = jet.abs() + swept.abs() + drive.body_flow_m3_s.abs() + wave.abs();
        let scale = self.stored_energy_j()
            + energy
            + aperture.pressure_work_j.abs()
            + mechanical.abs()
            + aperture.dissipated_energy_j.abs();
        if ![
            outgoing,
            aperture.bore_pressure_pa,
            jet,
            swept,
            wave,
            scale,
            flow_scale,
            aperture.balance_residual_j(),
            aperture.flow_residual_m3_s,
            aperture.dissipated_energy_j,
            mechanical,
        ]
        .iter()
        .all(|x| x.is_finite())
            || momentum.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "modal aperture refused an unresolved physical work or flow balance",
            ));
        }
        let (work_certificate, flow_certificate) =
            self.admit_balance(&aperture, scale, flow_scale, trial.evaluations)?;
        let volume = dt * (wave - drive.body_flow_m3_s);
        let plate_work = |column: usize| {
            dt * forces[column]
                * self
                    .spec
                    .modes
                    .iter()
                    .zip(&trial.velocity)
                    .map(|(m, v)| m.force_projection[column] * v)
                    .sum::<f64>()
        };
        let spatial_pressure = std::array::from_fn(|blade| {
            self.pressure_ports[blade].as_ref().map_or(
                pressure::MaterialPressureFrame::default(),
                |port| port.observe(pressures_pa[blade], &trial.velocity[..n], dt),
            )
        });
        if spatial_pressure.iter().any(|p| !p.flow_m3_s.is_finite() || !p.work_j.is_finite()) {
            return Err(invalid("material pressure work or conjugate flow overflowed"));
        }
        let port_transfers = [
            crate::port_adoption::PortTransfer::work(drive.upstream_pressure_pa * volume, true),
            crate::port_adoption::PortTransfer::work(-aperture.bore_pressure_pa * volume, true),
            crate::port_adoption::PortTransfer::work(plate_work(0), true),
            crate::port_adoption::PortTransfer::work(plate_work(1), true),
            crate::port_adoption::PortTransfer::work(spatial_pressure[0].work_j, true),
            crate::port_adoption::PortTransfer::work(spatial_pressure[1].work_j, true),
        ];
        Ok(ModalApertureFrame {
            aperture,
            modes: next,
            pressure_drop_pa: trial.dp,
            port_transfers,
            contact_force_n: contact.normal_force_n,
            momentum_residual_n: momentum,
            work_certificate,
            flow_certificate,
            spatial_pressure,
            memory: self
                .relaxation
                .as_ref()
                .map(|_| endpoint.expect("exact material endpoint").memory),
        })
    }

    fn admit_balance(
        &self,
        frame: &ApertureFrame,
        energy_scale: f64,
        flow_scale: f64,
        evaluations: usize,
    ) -> Result<(SolveCertificate, SolveCertificate), AcousticRealizeError> {
        let budget = SolveBudget {
            relative: self.spec.limits.relative_tolerance,
            absolute: 0.0,
            max_iterations: 0,
        };
        let work = budget.certificate(
            NewtonSite::ModalWork,
            frame
                .balance_residual_j()
                .abs()
                .max((-frame.dissipated_energy_j).max(0.0)),
            energy_scale.max(f64::MIN_POSITIVE),
            arithmetic_rounding_floor(4 * self.state.len() + 12, energy_scale),
            evaluations,
            0,
        );
        let flow = budget.certificate(
            NewtonSite::ModalFlow,
            frame.flow_residual_m3_s.abs(),
            flow_scale.max(f64::MIN_POSITIVE),
            arithmetic_rounding_floor(2 * self.state.len() + 8, flow_scale),
            evaluations,
            0,
        );
        Ok((
            admit_newton_iterate(work).map_err(AcousticRealizeError::AlgebraicSolve)?,
            admit_newton_iterate(flow).map_err(AcousticRealizeError::AlgebraicSolve)?,
        ))
    }
}

#[cfg(test)]
mod refusal_tests;

crate::port_adoption::impl_ported_module!(ModalAperture, ports);
