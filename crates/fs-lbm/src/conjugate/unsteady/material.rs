//! One native backward-Euler step on a retained, source-material momentum state.
//! Thermal feedback owns the clock and selects the accepted interval. This is
//! a checkpoint adapter around SIMPLEC, not another momentum discretization.

use super::{FlowStepRecord, UnsteadyConfig, admit_unsteady, converge_step};
use crate::conjugate::simple::{MomentumState, Solver, material_solver};
use crate::conjugate::{
    ChtError, FluidProperties, FvBoundary, FvFlow, SimpleConfig, TimeScheme, Turbulence,
    UnsteadyHeatingAudit, UnsteadyMechanicalStep, ViscousHeatControls, VoxelDomain, finite, poll,
};
use fs_exec::CancelGate;

#[derive(Clone)]
pub(in crate::conjugate) struct MaterialMomentum {
    pub time_s: f64,
    pub flow: FvFlow,
    pub record: Option<FlowStepRecord>,
    state: MomentumState,
    pub mechanical_step: Option<UnsteadyMechanicalStep>,
    heating: Option<MaterialHeating>,
}

#[derive(Clone)]
struct MaterialHeating {
    controls: ViscousHeatControls,
    audit: UnsteadyHeatingAudit,
    // Prepared by the accepted-duration hook, consumed by step only after
    // the original thermal owner and endpoint material checks succeed.
    pending: bool,
}

fn invalid(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "material_unsteady",
        reason: reason.into(),
    }
}

fn controls(dt: f64, config: &SimpleConfig) -> UnsteadyConfig {
    UnsteadyConfig {
        time_step_s: dt,
        steps: 1,
        scheme: TimeScheme::BackwardEuler,
        inner_iterations: config.max_iterations,
        inner_tolerance: config.tolerance,
        average_from_step: 1,
        probe: None,
    }
}

impl MaterialMomentum {
    pub fn at_rest(
        domain: &VoxelDomain,
        fluid: &FluidProperties,
        config: &SimpleConfig,
        time_s: f64,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        finite("material_unsteady.time_s", time_s)?;
        if config.turbulence != Turbulence::Laminar
            || config.fan.is_some()
            || !config.internal_fans.is_empty()
            || config
                .faces
                .iter()
                .any(|face| !matches!(face, FvBoundary::Wall { .. } | FvBoundary::Symmetry))
        {
            return Err(invalid(
                "retained material momentum requires laminar, closed wall/symmetry boundaries and no fans",
            ));
        }
        super::admit(domain, fluid, config)?;
        admit_unsteady(domain, config, &controls(1.0, config))?;
        let solver = Solver::new(domain, fluid, config);
        let state = solver.momentum_state();
        // Closed faces and native interior faces start at rest. No steady
        // solve and no manufactured velocity/pressure reconstruction occurs.
        let flow = solver.finish_projected(fluid, 0, (0.0, 0.0));
        poll(gate)?;
        Ok(Self {
            time_s,
            flow,
            record: None,
            state,
            mechanical_step: None,
            heating: None,
        })
    }

    pub fn has_heating(&self) -> bool {
        self.heating.is_some()
    }

    pub fn heating_audit(&self) -> Option<UnsteadyHeatingAudit> {
        self.heating.as_ref().map(|heating| heating.audit)
    }

    pub fn enable_heating(&mut self, controls: ViscousHeatControls) -> Result<(), ChtError> {
        if let Some(heating) = &mut self.heating {
            if heating.pending {
                return Err(invalid(
                    "change a transient actuator only at an accepted thermal endpoint",
                ));
            }
            heating.controls = controls;
        } else {
            if self.record.is_some() {
                return Err(invalid(
                    "select unsteady Stokes heating before the first momentum advance",
                ));
            }
            self.heating = Some(MaterialHeating {
                controls,
                audit: UnsteadyHeatingAudit::default(),
                pending: false,
            });
        }
        Ok(())
    }

    pub fn begin_advance(&mut self) -> Result<(), ChtError> {
        if let Some(heating) = &mut self.heating {
            if heating.pending {
                return Err(invalid("unfinished transient heat interval"));
            }
            let kinetic = self
                .mechanical_step
                .as_ref()
                .map_or(0.0, |step| step.final_kinetic_energy_j);
            heating.audit = UnsteadyHeatingAudit::from_kinetic_energy(kinetic);
        }
        Ok(())
    }

    pub fn accept_heating_step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
    ) -> Result<(), ChtError> {
        self.check_time(end_s)?;
        let heating = self
            .heating
            .as_mut()
            .ok_or_else(|| invalid("transient heating is not enabled"))?;
        let step = self
            .mechanical_step
            .as_ref()
            .ok_or_else(|| invalid("missing prepared transient heat"))?;
        if !heating.pending
            || step.start_s.to_bits() != start_s.to_bits()
            || step.duration_s.to_bits() != dt_s.to_bits()
            || step.end_s.to_bits() != end_s.to_bits()
        {
            return Err(invalid(
                "transient heat must accept its exact prepared thermal interval once",
            ));
        }
        heating.audit.accept(step, heating.controls)?;
        heating.pending = false;
        Ok(())
    }

    pub fn check_time(&self, time_s: f64) -> Result<(), ChtError> {
        if self.time_s.to_bits() != time_s.to_bits() {
            return Err(invalid(
                "retained momentum clock differs from this thermal state; resume its accepted coupled timeline",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &mut self,
        domain: &VoxelDomain,
        fluid: &FluidProperties,
        config: &SimpleConfig,
        viscosity_pa_s: &[f64],
        force_n_m3: Option<&[[f64; 3]]>,
        drag_pa_s_m2: Option<&[f64]>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        self.check_time(start_s)?;
        if self.heating.as_ref().is_some_and(|heating| heating.pending) {
            return Err(invalid(
                "transient momentum already prepared for this thermal interval",
            ));
        }
        // Deposition either adds dt or snaps to a stop whose exact difference
        // supplied dt. Admit precisely those two floating-point conventions.
        if !end_s.is_finite()
            || end_s <= start_s
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || ((start_s + dt_s).to_bits() != end_s.to_bits()
                && (end_s - start_s).to_bits() != dt_s.to_bits())
        {
            return Err(invalid(
                "momentum interval must equal the accepted thermal interval",
            ));
        }
        let unsteady = controls(dt_s, config);
        admit_unsteady(domain, config, &unsteady)?;
        let mut solver = material_solver(
            domain,
            fluid,
            viscosity_pa_s,
            force_n_m3,
            drag_pa_s_m2,
            config,
            gate,
        )?;
        if let Some(heating) = &self.heating {
            solver.enable_unsteady_stokes_heat(heating.controls)?;
        }
        solver.restore_momentum(&self.state)?;
        solver.start_time(dt_s, TimeScheme::BackwardEuler);
        let (inner_iterations, residuals) = converge_step(&mut solver, &unsteady, 1, gate)?;
        solver.project(gate)?;
        let mechanical_step = self
            .heating
            .as_ref()
            .map(|heating| solver.unsteady_mechanical_step(start_s, end_s, heating.controls, gate))
            .transpose()?;
        let mut kinetic_energy_j = 0.0;
        for c in 0..domain.cell_count() {
            if c % 256 == 0 {
                poll(gate)?;
            }
            if domain.is_fluid(c) {
                let v = solver.cell_velocity(c);
                kinetic_energy_j += 0.5
                    * fluid.density_kg_m3
                    * domain.volume(c)
                    * v.iter().map(|component| component * component).sum::<f64>();
            }
        }
        finite("material_unsteady.kinetic_energy_j", kinetic_energy_j)?;
        let state = solver.momentum_state();
        let flow = solver.finish_projected(fluid, inner_iterations, residuals);
        let record = FlowStepRecord {
            time_s: end_s,
            inner_iterations,
            mass_residual: residuals.0,
            momentum_residual: residuals.1,
            probe_velocity_m_s: [0.0; 3],
            kinetic_energy_j,
        };
        poll(gate)?;
        self.state = state;
        self.flow = flow;
        self.record = Some(record);
        self.mechanical_step = mechanical_step;
        if let Some(heating) = &mut self.heating {
            heating.pending = true;
        }
        self.time_s = end_s;
        Ok(())
    }
}
