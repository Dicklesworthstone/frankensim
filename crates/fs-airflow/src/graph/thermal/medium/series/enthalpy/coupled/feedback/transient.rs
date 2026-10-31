//! Finite pipe inventory and retained momentum for a closed Boussinesq loop.
//!
//! The original caloric cells store heat even at zero flow. The original PHS
//! inertance stores momentum; a reversal never replaces it with a steady root.
//! Density varies ONLY in gravity. Reference mass, volume, continuity and inertia
//! use rho0; viscosity is resolved from the same bound source at every trial.
//!
//! This is a low-Eckert, single-phase approximation. Thermal and mechanical
//! accounts are separate: neither pressure work nor viscous/gravity work is
//! fabricated as a heat source. Non-cancelling omitted-work bounds are explicit.
//! No compressible total-energy, pipe-wall film, boiling or turbulent closure is
//! inferred. Native thermal links/boundaries remain active while flow is zero.

/// Source-backed fixed-composition gas circulation on the same finite inventory.
pub mod gas;
mod step;
pub use gas::{FiniteGasLoop, GasLoopLimits};

use super::super::super::{LiquidEnthalpyError, admit_liquid_caloric_source};
use crate::graph::thermal::medium::series::hydraulic::{
    MaterialFlowResistance, ResistanceError, ResistanceEvaluation,
    pressure::{HydrostaticPressure, SeriesGravity, SeriesPressureError},
};
use fs_conduction::transient::deposition::{
    CellGeometry, CellId, DepositionError, DepositionSimulation, EnergyAudit, MaterialId,
    ThermalLaw,
};
use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::MaterialStatePointError;
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, StepRecord};
use fs_qty::{MassFlowRate, Temperature};
use fs_solver::{NewtonError, NewtonKrylovConfig, NewtonReport};
use std::collections::BTreeSet;
use std::convert::Infallible;

/// Original finite fluid cell and its authored linear/quadratic loss geometry.
/// Order is the positive circulation direction, independent of actual flow sign.
#[derive(Clone, Copy, Debug)]
pub struct PipeCell {
    pub cell: CellId,
    pub resistance: MaterialFlowResistance,
}

/// Explicit spatial/model admission, independent of nonlinear tolerances.
#[derive(Clone, Copy, Debug)]
pub struct FiniteLoopLimits {
    /// Includes any finite wall cells sharing the native thermal residual.
    pub maximum_cells: usize,
    /// Maximum |rho(T)/rho0 - 1| at EVERY constitutive trial, less than one.
    pub maximum_relative_density_departure: f64,
}

/// Controls on one joint solve. Algebraic convergence does not certify time or
/// upwind-cell resolution; refine both against the physical response.
#[derive(Clone, Copy, Debug)]
pub struct FiniteLoopControls {
    pub solver: NewtonKrylovConfig,
    pub maximum_newton_iterations: usize,
    /// Counts every complete residual proposal, including numerical differences
    /// and final acceptance evaluation. No source error is a missing flow root.
    pub maximum_residual_evaluations: usize,
    /// Positive coordinate scales; the solver acts on dimensionless residuals.
    pub energy_scale_j: f64,
    pub flow_scale_m3_s: f64,
    pub finite_difference_step: f64,
    pub cell_tolerance_j: f64,
    pub total_thermal_tolerance_j: f64,
    pub momentum_tolerance_pa_s: f64,
    pub mechanical_energy_tolerance_j: f64,
    /// Hard low-Eckert caps using the largest absolute old/new/average flow,
    /// with the same held endpoint constitutive states. This bounds the admitted
    /// discrete interval approximation, not unresolved continuous-time extrema.
    /// Loss and absolute LOCAL anomalous-gravity work cannot cancel each other.
    pub maximum_neglected_power_w: f64,
    pub maximum_neglected_energy_j: f64,
    /// Maximum dt*R/L, in (0, 2]. Prevents passive midpoint drag from producing
    /// an unforced sign flip; this is not a temporal-error certificate.
    pub maximum_drag_step: f64,
}

/// Original refusal, numerical failure or explicit approximation limit.
#[derive(Debug)]
pub enum FiniteLoopError<E = Infallible> {
    Invalid(&'static str),
    Caloric(LiquidEnthalpyError),
    Material(MaterialStatePointError),
    GasCaloric(fs_material::thermochemical_gas::energy::GasEnergyError),
    GasMach {
        section: usize,
        value: f64,
        maximum: f64,
    },
    Transport(MaterialError),
    Thermal(DepositionError),
    Resistance {
        section: usize,
        source: ResistanceError,
    },
    Gravity(SeriesPressureError),
    Momentum(PhsError),
    Solver(NewtonError),
    NotConverged(NewtonReport),
    Heat(E),
    DensityDeparture {
        section: usize,
        relative: f64,
        maximum: f64,
    },
    NeglectedWork {
        power_w: f64,
        energy_j: f64,
    },
    DragStep {
        value: f64,
        maximum: f64,
    },
    WorkBudget {
        performed: usize,
        maximum: usize,
    },
    MomentumResidual {
        impulse_pa_s: f64,
        energy_j: f64,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for FiniteLoopError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite liquid loop: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for FiniteLoopError<E> {}
impl<E> From<DepositionError> for FiniteLoopError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// Source and original loss evaluations at the ACTUAL native step-average flow.
#[derive(Debug)]
pub struct FiniteLoopHydraulics<F = ResolvedLiquidState> {
    pub sources: Vec<F>,
    pub losses: Vec<ResistanceEvaluation>,
    /// Original loss/Reynolds admission at the nonnegative speed envelope.
    /// These bound omitted work; they are not native discrete dissipation.
    pub envelope_losses: Vec<ResistanceEvaluation>,
    pub hydrostatic: HydrostaticPressure,
    pub momentum: StepRecord,
    pub initial_volume_flow_m3_s: f64,
    pub final_volume_flow_m3_s: f64,
    /// Used by every advective link. During reversal its sign can differ from
    /// the final flow sign; endpoint flow is never substituted for this value.
    pub average_volume_flow_m3_s: f64,
    pub driving_pressure_pa: f64,
    pub buoyancy_work_j: f64,
    /// Original loss at the signed native average flow, integrated over dt.
    pub hydraulic_loss_j: f64,
    /// max(|Q_old|, |Q_new|, |Q_average|), with endpoint-held source states.
    pub neglected_work_speed_envelope_m3_s: f64,
    pub neglected_hydraulic_energy_bound_j: f64,
    pub neglected_buoyancy_exchange_bound_j: f64,
    pub neglected_energy_bound_j: f64,
    /// Physical original-law residual, independently of the nonlinear stopping norm.
    pub momentum_residual_pa_s: f64,
    /// Delta kinetic energy + original loss - external pressure work - buoyancy work.
    pub mechanical_balance_residual_j: f64,
}

/// A successful joint publication, with the exact supplied heat-source payload.
#[derive(Debug)]
pub struct FiniteLoopStep<T, F = ResolvedLiquidState> {
    pub hydraulics: FiniteLoopHydraulics<F>,
    pub heat: T,
    pub thermal: EnergyAudit,
    pub maximum_cell_residual_j: f64,
    pub residual_evaluations: usize,
    pub nonlinear: NewtonReport,
}

#[derive(Clone, Debug)]
struct BoundCell {
    section: PipeCell,
    material: MaterialId,
    mass_kg: f64,
    geometry: CellGeometry,
    specific_gravity_work_j_kg: f64,
}

/// One accepted hydraulic momentum and the original shared thermal timeline.
/// The thermal owner supplies every time interval. The retained accepted time
/// only prevents attaching old momentum to an independently advanced timeline.
pub struct FiniteLiquidLoop<'a, 'card> {
    law: &'a LiquidTemperatureLaw<'card>,
    source: &'a ResolvedHeatCapacityCurve,
    core: FiniteLoopCore<'a>,
}

// The geometry, one momentum and native thermal transaction are shared; source
// admission stays with the concrete liquid/gas adapters. This is not a second
// material interface or a solver.
struct FiniteLoopCore<'a> {
    caloric: &'a dyn ThermalLaw,
    cells: Vec<BoundCell>,
    gravity: SeriesGravity,
    limits: FiniteLoopLimits,
    reference_density: f64,
    inertance: f64,
    momentum_pa_s: f64,
    accepted_time_s: f64,
}

impl<'a, 'card> FiniteLiquidLoop<'a, 'card> {
    /// Bind actual cells to the original liquid source and caloric object.
    /// Geometry, reference inventory and source admission are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cx: &Cx<'_>,
        law: &'a LiquidTemperatureLaw<'card>,
        source: &'a ResolvedHeatCapacityCurve,
        thermal: &DepositionSimulation<'_>,
        sections: &[PipeCell],
        gravity: SeriesGravity,
        initial_volume_flow_m3_s: f64,
        limits: FiniteLoopLimits,
    ) -> Result<Self, FiniteLoopError> {
        poll(cx)?;
        admit_liquid_caloric_source(cx, law.reference_state(), source, law.temperature_axis())
            .map_err(FiniteLoopError::Caloric)?;
        let reference = law
            .reference_state()
            .transport_state()
            .map_err(FiniteLoopError::Transport)?;
        let core = FiniteLoopCore::try_new(
            cx,
            source.curve(),
            reference,
            thermal,
            sections,
            gravity,
            initial_volume_flow_m3_s,
            limits,
        )?;
        let (_, transports) = core.resolve(cx, thermal, &|_, temperature| {
            law.at_temperature(temperature)
                .map_err(FiniteLoopError::<Infallible>::Material)
        })?;
        core.losses::<Infallible>(cx, &transports, initial_volume_flow_m3_s)?;
        Ok(Self { law, source, core })
    }
    pub const fn source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.source
    }
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.core.reference_density
    }
    pub const fn inertance_pa_s2_m3(&self) -> f64 {
        self.core.inertance
    }
    pub const fn accepted_time_s(&self) -> f64 {
        self.core.accepted_time_s
    }
    pub fn volume_flow_m3_s(&self) -> f64 {
        self.core.volume_flow_m3_s()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.core.kinetic_energy_j()
    }
}

impl<'a> FiniteLoopCore<'a> {
    #[allow(clippy::too_many_arguments)]
    fn try_new(
        cx: &Cx<'_>,
        caloric: &'a dyn ThermalLaw,
        reference: FluidTransportState,
        thermal: &DepositionSimulation<'_>,
        sections: &[PipeCell],
        gravity: SeriesGravity,
        initial_volume_flow_m3_s: f64,
        limits: FiniteLoopLimits,
    ) -> Result<Self, FiniteLoopError> {
        poll(cx)?;
        if !(2..=4096).contains(&limits.maximum_cells)
            || sections.len() < 2
            || sections.len() > limits.maximum_cells
            || thermal.cells().len() > limits.maximum_cells
            || !nonnegative(limits.maximum_relative_density_departure)
            || limits.maximum_relative_density_departure >= 1.0
            || !initial_volume_flow_m3_s.is_finite()
            || gravity.section_count() != sections.len()
            || gravity.vertices_m().first() != gravity.vertices_m().last()
        {
            return Err(FiniteLoopError::Invalid(
                "bounded closed finite-cell geometry and Boussinesq limits required",
            ));
        }
        let rho0 = reference.density_kg_m3();
        let mut seen = BTreeSet::new();
        let mut cells = Vec::with_capacity(sections.len());
        let mut inertance = 0.0;
        for (i, section) in sections.iter().enumerate() {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&section.cell)
                .ok_or(FiniteLoopError::Invalid("pipe cell is not active"))?;
            if !seen.insert(section.cell)
                || !thermal.uses_material_law(cell.material_id, caloric)
                || cell.geometry.start_m != gravity.vertices_m()[i]
                || cell.geometry.end_m != gravity.vertices_m()[i + 1]
            {
                return Err(FiniteLoopError::Invalid(
                    "unique cells, original caloric law and exact gravity endpoints required",
                ));
            }
            let delta: [f64; 3] =
                core::array::from_fn(|a| cell.geometry.end_m[a] - cell.geometry.start_m[a]);
            let length = fs_math::det::sqrt(delta.iter().map(|v| v * v).sum());
            let area = section.resistance.flow_area().value();
            if !positive(length)
                || !near(cell.geometry.volume_m3, area * length)
                || !near(cell.mass_kg, rho0 * cell.geometry.volume_m3)
            {
                return Err(FiniteLoopError::Invalid(
                    "pipe inventory must satisfy V=A*length and M=rho0*V",
                ));
            }
            inertance += rho0 * length / area;
            let work = delta
                .iter()
                .zip(gravity.gravity_m_s2())
                .map(|(dx, g)| dx * g)
                .sum::<f64>();
            if !work.is_finite() || !positive(inertance) {
                return Err(FiniteLoopError::Invalid(
                    "pipe inertance or gravity work is unrepresentable",
                ));
            }
            section
                .resistance
                .evaluate(&reference, MassFlowRate::new(0.0))
                .map_err(|source| FiniteLoopError::Resistance { section: i, source })?;
            cells.push(BoundCell {
                section: *section,
                material: cell.material_id,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                specific_gravity_work_j_kg: work,
            });
        }
        let momentum = inertance * initial_volume_flow_m3_s;
        if !momentum.is_finite()
            || !positive(1.0 / inertance)
            || !(0.5 * momentum * initial_volume_flow_m3_s).is_finite()
            || (initial_volume_flow_m3_s != 0.0 && momentum == 0.0)
        {
            return Err(FiniteLoopError::Invalid(
                "initial pipe momentum is unrepresentable",
            ));
        }
        let owner = Self {
            caloric,
            cells,
            gravity,
            limits,
            reference_density: rho0,
            inertance,
            momentum_pa_s: momentum,
            accepted_time_s: thermal.time_s(),
        };
        Ok(owner)
    }

    fn volume_flow_m3_s(&self) -> f64 {
        self.momentum_pa_s / self.inertance
    }
    fn kinetic_energy_j(&self) -> f64 {
        0.5 * self.momentum_pa_s * self.volume_flow_m3_s()
    }

    fn resolve<E, F: FluidTransport>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        resolve: &impl Fn(CellId, Temperature) -> Result<F, FiniteLoopError<E>>,
    ) -> Result<(Vec<F>, Vec<FluidTransportState>), FiniteLoopError<E>> {
        if thermal.cells().len() > self.limits.maximum_cells {
            return Err(FiniteLoopError::Invalid(
                "finite loop thermal cell budget exceeded",
            ));
        }
        let mut sources = Vec::with_capacity(self.cells.len());
        let mut transports = Vec::with_capacity(self.cells.len());
        for (i, bound) in self.cells.iter().enumerate() {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&bound.section.cell)
                .ok_or(FiniteLoopError::Invalid("bound pipe cell disappeared"))?;
            if cell.material_id != bound.material
                || cell.mass_kg != bound.mass_kg
                || cell.geometry != bound.geometry
                || !thermal.uses_material_law(cell.material_id, self.caloric)
            {
                return Err(FiniteLoopError::Invalid(
                    "pipe geometry, mass or original caloric source changed",
                ));
            }
            let temperature = thermal.temperature_k(bound.section.cell)?;
            let source = resolve(bound.section.cell, Temperature::new(temperature))?;
            let transport = source
                .transport_state()
                .map_err(FiniteLoopError::Transport)?;
            let relative = (transport.density_kg_m3() / self.reference_density - 1.0).abs();
            if !relative.is_finite() || relative > self.limits.maximum_relative_density_departure {
                return Err(FiniteLoopError::DensityDeparture {
                    section: i,
                    relative,
                    maximum: self.limits.maximum_relative_density_departure,
                });
            }
            sources.push(source);
            transports.push(transport);
        }
        poll(cx)?;
        Ok((sources, transports))
    }

    fn losses<E>(
        &self,
        cx: &Cx<'_>,
        sources: &[FluidTransportState],
        flow: f64,
    ) -> Result<Vec<ResistanceEvaluation>, FiniteLoopError<E>> {
        let mass = self.reference_density * flow;
        if !mass.is_finite() || (flow != 0.0 && mass == 0.0) {
            return Err(FiniteLoopError::Invalid(
                "reference mass flow is unrepresentable",
            ));
        }
        self.cells
            .iter()
            .zip(sources)
            .enumerate()
            .map(|(i, (cell, source))| {
                poll(cx)?;
                // Explicit Boussinesq numerical view, beside the untouched source:
                // actual source rho supplies gravity; only rho0 supplies continuity.
                let reference = FluidTransportState::try_new(
                    self.reference_density,
                    source.specific_heat_j_kg_k(),
                    source.thermal_conductivity_w_m_k(),
                    source.dynamic_viscosity_pa_s(),
                )
                .map_err(FiniteLoopError::Transport)?;
                cell.section
                    .resistance
                    .evaluate(&reference, MassFlowRate::new(mass))
                    .map_err(|source| FiniteLoopError::Resistance { section: i, source })
            })
            .collect()
    }

    fn hydraulics<E, F: FluidTransport>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        resolve: &impl Fn(CellId, Temperature) -> Result<F, FiniteLoopError<E>>,
        candidate_average_flow: f64,
        drive: f64,
        dt: f64,
        controls: FiniteLoopControls,
    ) -> Result<FiniteLoopHydraulics<F>, FiniteLoopError<E>> {
        let (sources, transports) = self.resolve(cx, thermal, resolve)?;
        let trial_losses = self.losses(cx, &transports, candidate_average_flow)?;
        let resistance: f64 = trial_losses
            .iter()
            .map(|loss| {
                if candidate_average_flow == 0.0 {
                    loss.dpressure_dmass_flow * self.reference_density
                } else {
                    loss.pressure_drop_pa / candidate_average_flow
                }
            })
            .sum();
        if !nonnegative(resistance) {
            return Err(FiniteLoopError::Invalid(
                "nonnegative finite effective pipe resistance required",
            ));
        }
        let drag_step = dt * (resistance / self.inertance);
        if !drag_step.is_finite() || drag_step > controls.maximum_drag_step {
            return Err(FiniteLoopError::DragStep {
                value: drag_step,
                maximum: controls.maximum_drag_step,
            });
        }
        let hydrostatic = self
            .gravity
            .evaluate_densities(
                transports.iter().map(FluidTransportState::density_kg_m3),
                0.0,
                || cx.checkpoint().map_err(|_| SeriesPressureError::Cancelled),
            )
            .map_err(FiniteLoopError::Gravity)?;
        let pressure = drive + hydrostatic.pressure_gain_pa();
        if !pressure.is_finite() {
            return Err(FiniteLoopError::Invalid(
                "total pressure drive is unrepresentable",
            ));
        }
        let storage = QuadraticStorage::new(vec![1.0 / self.inertance], 1)
            .map_err(FiniteLoopError::Momentum)?;
        let system = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![resistance],
            vec![1.0],
            Box::new(storage),
        )
        .map_err(FiniteLoopError::Momentum)?;
        poll(cx)?;
        let momentum = fs_phs::step(&system, &[self.momentum_pa_s], &[pressure], dt)
            .map_err(FiniteLoopError::Momentum)?;
        poll(cx)?;
        let average = momentum.y[0];
        let final_flow = momentum.x[0] / self.inertance;
        // Re-evaluate the ORIGINAL nonlinear loss at the actual native flow.
        // Acceptance below gates its own impulse/power residual, so a lagged
        // candidate coefficient cannot masquerade as a closed momentum solve.
        let losses = self.losses(cx, &transports, average)?;
        let loss_pressure: f64 = losses.iter().map(|loss| loss.pressure_drop_pa).sum();
        let loss_power: f64 = losses.iter().map(|loss| loss.dissipated_power_w).sum();
        let actual_resistance = if average == 0.0 {
            losses
                .iter()
                .map(|loss| loss.dpressure_dmass_flow * self.reference_density)
                .sum()
        } else {
            loss_pressure / average
        };
        let actual_drag_step = dt * (actual_resistance / self.inertance);
        if !nonnegative(actual_resistance)
            || !actual_drag_step.is_finite()
            || actual_drag_step > controls.maximum_drag_step
        {
            return Err(FiniteLoopError::DragStep {
                value: actual_drag_step,
                maximum: controls.maximum_drag_step,
            });
        }
        // Signed average flow remains authoritative for native transport/work.
        // It can vanish during reversal despite nonzero speed throughout most
        // of the interval. Bound omissions and admit Reynolds at the largest
        // endpoint/average speed under the SAME held constitutive states.
        let speed_envelope = self
            .volume_flow_m3_s()
            .abs()
            .max(final_flow.abs())
            .max(average.abs());
        let envelope_losses = self.losses(cx, &transports, speed_envelope)?;
        let envelope_loss_power: f64 = envelope_losses
            .iter()
            .map(|loss| loss.dissipated_power_w)
            .sum();
        let buoyancy_power = hydrostatic.pressure_gain_pa() * average;
        let local_buoyancy_bound: f64 = self
            .cells
            .iter()
            .zip(&transports)
            .map(|(cell, source)| {
                ((source.density_kg_m3() - self.reference_density)
                    * cell.specific_gravity_work_j_kg
                    * speed_envelope)
                    .abs()
            })
            .sum();
        let neglected_power = envelope_loss_power + local_buoyancy_bound;
        let neglected_energy = dt * neglected_power;
        let impulse = (momentum.x[0] - self.momentum_pa_s) - dt * (pressure - loss_pressure);
        let loss_energy = dt * loss_power;
        let buoyancy_work = dt * buoyancy_power;
        let energy_residual = momentum.delta_h + loss_energy - dt * drive * average - buoyancy_work;
        if [
            average,
            final_flow,
            speed_envelope,
            loss_pressure,
            loss_power,
            envelope_loss_power,
            buoyancy_power,
            local_buoyancy_bound,
            neglected_power,
            neglected_energy,
            impulse,
            energy_residual,
            loss_energy,
            buoyancy_work,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(FiniteLoopError::Invalid(
                "finite-loop work or momentum is unrepresentable",
            ));
        }
        if neglected_power > controls.maximum_neglected_power_w
            || neglected_energy > controls.maximum_neglected_energy_j
        {
            return Err(FiniteLoopError::NeglectedWork {
                power_w: neglected_power,
                energy_j: neglected_energy,
            });
        }
        Ok(FiniteLoopHydraulics {
            sources,
            losses,
            envelope_losses,
            hydrostatic,
            momentum,
            initial_volume_flow_m3_s: self.volume_flow_m3_s(),
            final_volume_flow_m3_s: final_flow,
            average_volume_flow_m3_s: average,
            driving_pressure_pa: drive,
            buoyancy_work_j: buoyancy_work,
            hydraulic_loss_j: loss_energy,
            neglected_work_speed_envelope_m3_s: speed_envelope,
            neglected_hydraulic_energy_bound_j: dt * envelope_loss_power,
            neglected_buoyancy_exchange_bound_j: dt * local_buoyancy_bound,
            neglected_energy_bound_j: neglected_energy,
            momentum_residual_pa_s: impulse,
            mechanical_balance_residual_j: energy_residual,
        })
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), FiniteLoopError<E>> {
    cx.checkpoint().map_err(|_| FiniteLoopError::Cancelled)
}
fn thermal_poll(cx: &Cx<'_>) -> Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn nonnegative(v: f64) -> bool {
    v.is_finite() && v >= 0.0
}
fn near(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}
