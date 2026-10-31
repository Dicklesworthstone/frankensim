//! Buoyancy-selected circulation on the original finite-cell heat trajectory.
//!
//! Fixed inventories use reference density while their source-resolved bulk
//! temperatures set the hydrostatic contrast. The original implicit caloric
//! residual owns every advected joule and reciprocal wall transfer. This is a
//! low-Eckert approximation: pressure/expansion work and viscous heat are not
//! converted into thermal storage, and their explicit W/J bounds must hold.

use super::*;
use crate::graph::thermal::medium::series::hydraulic::pressure::{
    BoussinesqDensity, SeriesGravity,
};

/// A closed path through distinct original coolant cells. Vertices follow
/// section order and must return to the exact starting point. Section density
/// contrasts follow the bulk state, never the wall/film temperature.
#[derive(Clone, Copy, Debug)]
pub struct GravityCoolantLoop<'a> {
    pub sections: &'a [PressureCoolantSection],
    pub links: &'a [PressureConvectiveLink],
    pub gravity: &'a SeriesGravity,
    pub density: BoussinesqDensity,
    /// Additional signed differential pressure [Pa], including zero.
    pub driving_pressure_pa: f64,
    pub pressure_budget: SeriesPressureBudget,
    pub maximum_links: usize,
    /// Bound on the sum of absolute imposed-pressure work, buoyancy work and
    /// viscous loss [W]. Absolute terms prevent opposing drives from hiding
    /// a large omitted mechanical/thermal conversion behind a small net loss.
    pub maximum_omitted_conversion_w: f64,
    /// Bound on that same sum integrated over this accepted interval [J].
    pub maximum_omitted_conversion_j: f64,
}

/// Mechanical quantities deliberately excluded from the caloric equation.
/// These are physical omission estimates, not thermal residual tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoolantGravityWork {
    /// Signed external differential-pressure power `dp_imposed * Q0` [W].
    pub imposed_pressure_w: f64,
    /// Signed density-contrast power `head * Q0` [W]. This is NOT heater work
    /// and is distinct from conservative `mdot*g.dot(last-first)` (zero here).
    pub buoyancy_w: f64,
    /// Actual positive hydraulic loss, not deposited a second time as heat [W].
    pub dissipation_w: f64,
    /// Independently evaluated imposed + buoyancy - loss [W].
    pub balance_residual_w: f64,
    pub omitted_absolute_conversion_w: f64,
    pub omitted_absolute_conversion_j: f64,
}

/// Source owners and mechanical omission at the EXACT accepted caloric trial.
#[derive(Debug)]
pub struct GravityRelaxationConvection<F, L = F> {
    pub hydraulics: SeriesPressureFlow<L>,
    pub convection: FiniteRelaxationConvection<F, L>,
    /// Local losses excluded under the declared low-Eckert work budget [W].
    pub omitted_hydraulic_heat_w: BTreeMap<CellId, f64>,
    pub work: CoolantGravityWork,
}

impl<R> SourceConvectionHeatNetwork<'_, '_, '_, '_, R> {
    /// Evolve closed finite coolant inventories with signed natural circulation.
    ///
    /// Every native caloric trial re-resolves the original bulk rho/mu, solves
    /// the original resistance inverse, updates wall-transfer speed, and feeds
    /// the SAME flow into conservative upwind enthalpy exchange. Reversal and
    /// hydrostatic rest need no separate timeline. Geometry, stored mass and
    /// reference volume stay fixed; actual density departures must remain in
    /// the caller's Boussinesq bound and each cell must have reference mass/V.
    ///
    /// Thermal storage, circuit energy and material memory remain with the
    /// existing transactional owner. No gravity power is mislabeled as external
    /// heat. Both pressure/expansion conversion and viscous heating are omitted
    /// under explicit absolute W/J limits; exact total thermodynamic energy is
    /// not claimed. The native source/heat/cancellation gates still apply.
    ///
    /// Supplied links retain their original correlation support. In particular,
    /// a forced-convection card may refuse a stagnant state; do not invent a
    /// minimum velocity. Empty links admit an adiabatic loop, including rest.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_gravity_driven_relaxation_with_links<
        F: ThermalFluidTransport,
        L: ThermalFluidTransport,
        M,
        E,
    >(
        &self,
        cx: &Cx<'_>,
        owner: &mut FrameRelaxationThermalSystem<'_, '_, '_, '_, M>,
        thermal: &mut DepositionSimulation<'_>,
        external: &[f64],
        electrical: Option<RelaxationElectricalInput<'_>>,
        controls: FrameThermalSolve,
        energy_guess_j: &[f64],
        gate: &CancelGate,
        path: GravityCoolantLoop<'_>,
        resolve_coolant: impl Fn(
            &PressureCoolantSection,
            &DepositionSimulation<'_>,
            Temperature,
        ) -> Result<L, E>,
    ) -> Result<
        RelaxationHeatReport<GravityRelaxationConvection<F, L>>,
        RelaxationTransportError<PressureCoolantError<E>>,
    >
    where
        R: Fn(&SourceConvectiveCell, &DepositionSimulation<'_>, Temperature) -> Result<F, E>,
    {
        checkpoint(cx, gate).map_err(RelaxationTransportError::Source)?;
        let vertices = path.gravity.vertices_m();
        if path.gravity.section_count() != path.sections.len()
            || vertices.first() != vertices.last()
            || !path.driving_pressure_pa.is_finite()
            || [
                path.maximum_omitted_conversion_w,
                path.maximum_omitted_conversion_j,
            ]
            .into_iter()
            .any(|value| !value.is_finite() || value < 0.0)
        {
            return Err(RelaxationTransportError::Source(
                PressureCoolantError::Invalid(
                    "closed gravity geometry, finite drive and nonnegative omitted-work W/J bounds required",
                ),
            ));
        }
        let dt = owner.mechanical().sample_period_s();
        self.step_driven_relaxation(cx, owner, thermal, external, electrical, controls,
            energy_guess_j, gate, path.sections, path.links, path.pressure_budget.maximum_sections,
            path.maximum_links, false, resolve_coolant, |sections| {
                let flow = SeriesPressureFlow::solve_with_boussinesq_gravity(sections,
                    path.driving_pressure_pa, path.gravity.clone(), path.density,
                    path.pressure_budget, gate).map_err(PressureCoolantError::Hydraulic)?;
                let q = flow.evaluations()[0].volume_flow_m3_s;
                let imposed_pressure_w = finite(path.driving_pressure_pa * q)?;
                let buoyancy_w = finite(flow.hydrostatic().expect("gravity solve retains head")
                    .pressure_gain_pa() * q)?;
                let dissipation_w = flow.dissipated_power_w();
                let omitted_absolute_conversion_w = finite(imposed_pressure_w.abs()
                    + buoyancy_w.abs() + dissipation_w)?;
                let omitted_absolute_conversion_j = finite(dt * omitted_absolute_conversion_w)?;
                if omitted_absolute_conversion_w > path.maximum_omitted_conversion_w
                    || omitted_absolute_conversion_j > path.maximum_omitted_conversion_j {
                    return Err(PressureCoolantError::Invalid(
                        "gravity-loop mechanical/thermal conversion exceeds the declared W/J omission bound"));
                }
                let work = CoolantGravityWork { imposed_pressure_w, buoyancy_w, dissipation_w,
                    balance_residual_w: finite(imposed_pressure_w + buoyancy_w - dissipation_w)?,
                    omitted_absolute_conversion_w, omitted_absolute_conversion_j };
                checkpoint(cx, gate)?;
                Ok((flow, work))
            }, |hydraulics, convection, omitted_hydraulic_heat_w, work| {
                GravityRelaxationConvection { hydraulics, convection, omitted_hydraulic_heat_w, work }
            })
    }
}
