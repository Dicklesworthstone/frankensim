//! Implicit dynamic hydraulic loads on the original BDF Stokes operator.

use super::{
    StokesAssembly, StokesBoundaryError, StokesResistancePort, StokesSolution, invalid,
    velocity_port_flow,
};
use crate::ns::boundary::{StokesBoundary, bdf::StokesBdfStep};
use crate::ns::{FluxParams, FluxSystem};

#[path = "windkessel/model.rs"]
mod model;
#[path = "windkessel/nonlinear.rs"]
mod nonlinear;
use model::ReducedWindkesselPort;
pub use model::{StokesWindkesselBalance, StokesWindkesselPort};

/// One common time level for fluid momentum and external compliance storage.
/// Port order and fixed parameters/reference must match the new step.
#[derive(Clone, Copy)]
pub struct StokesWindkesselHistory<'a> {
    /// Original BDM velocity moments, excluding pressure coefficients [m/s].
    /// All flow histories are derived from these moments at the actual ports;
    /// a detached independently supplied inertance-flow history is not accepted.
    pub velocity: &'a [f64],
    /// Compliance pressures relative to each reservoir [m²/s²], in port order.
    /// These are the preceding accepted balance's compliance-pressure values.
    pub compliance_pressures_m2_s2: &'a [f64],
}

/// Accepted native field plus the dynamic state and power of every external load.
pub struct StokesWindkesselSolution<S = StokesSolution> {
    /// Original continuum field, pressure datum and discrete boundary reactions.
    pub flow: S,
    /// Resolved-fluid kinetic energy divided by reference density [m⁴/s²],
    /// evaluated using the same retained mass metric as the momentum step.
    pub kinetic_energy_m4_s2: f64,
    /// Added uniform port pressure increments [m²/s²].
    pub pressure_offsets_m2_s2: Vec<f64>,
    /// Actual signed outward integral flow per unit depth [m²/s].
    pub outward_volume_rates_m2_s: Vec<f64>,
    /// Accepted compliance coordinates, storage, resistor losses and BDF defect.
    pub balances: Vec<StokesWindkesselBalance>,
    /// Largest absolute original pressure/inertance-law residual [m²/s²].
    pub maximum_pressure_residual_m2_s2: f64,
    /// Largest absolute original compliance-continuity residual [m²/s].
    pub maximum_flow_residual_m2_s: f64,
}

impl FluxSystem<'_> {
    /// Advance laminar Stokes and passive RCR loads at one common implicit endpoint.
    ///
    /// For outward Q and compliance pressure x relative to the held reservoir:
    /// `p=p_res+Rp*Q+L*D_t(Q)+x`, `C*D_t(x)=Q-x/Rd`. The same original BDF
    /// derivative advances fluid momentum, external inertance and compliance.
    /// `older_history=None` uses BDF1; otherwise its interval belongs to BOTH
    /// older fluid and capacitor states, with the original step-ratio cap of 2.
    /// Every previous/older Q is integrated from its native velocity moments.
    ///
    /// Eliminating x and history gives one nonnegative endpoint resistance per
    /// port. The existing native pressure responses solve all ports together,
    /// retaining resolved profiles, cross-port coupling and original tractions.
    /// Every traction port may be controlled; reservoir laws select pressure.
    /// Parameters and the compliance reference are fixed through the histories;
    /// reservoir increments are endpoint data held over the current interval.
    /// Pressure power refers to the added -p*n load; original traction work is
    /// separate. This does not transport density, heat or an external interface.
    ///
    /// Both original dynamic laws and native continuum residuals are checked.
    /// External stored energy and resistor losses are separate from fluid heat;
    /// the signed BDF energy defect is reported, not reclassified as dissipation.
    /// No arbitrary-step passivity, nonlinear rheology, cavitation or distributed
    /// wave propagation claim follows from this lumped linear load.
    ///
    /// # Errors
    /// Refuses invalid/incomplete histories, intervals, passive parameters,
    /// normalization, port ownership, singular responses or failed equations.
    /// Histories remain immutable; no partially updated load state is published.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_stokes_mixed_variable_step_windkessel_ports(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[StokesWindkesselPort],
        previous: StokesWindkesselHistory<'_>,
        older_history: Option<(StokesWindkesselHistory<'_>, f64)>,
        timestep_s: f64,
        absolute_pressure_tolerance_m2_s2: f64,
        absolute_flow_tolerance_m2_s: f64,
    ) -> Result<StokesWindkesselSolution, StokesBoundaryError> {
        let assembly = self.assemble_stokes_mixed_bdf2_variable_step(
            params,
            force,
            boundaries,
            previous.velocity,
            older_history.map(|(history, dt)| (history.velocity, dt)),
            timestep_s,
        )?;
        solve_windkessel_assembly(
            assembly,
            ports,
            loads,
            previous,
            older_history,
            timestep_s,
            absolute_pressure_tolerance_m2_s2,
            absolute_flow_tolerance_m2_s,
        )
    }
}

/// Only the common native and material BDF constructors call this entry; an
/// arbitrary prepared steady/frozen operator cannot be paired with a new dt.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_windkessel_assembly(
    assembly: StokesAssembly<'_, '_>,
    ports: &[&[usize]],
    loads: &[StokesWindkesselPort],
    previous: StokesWindkesselHistory<'_>,
    older_history: Option<(StokesWindkesselHistory<'_>, f64)>,
    timestep_s: f64,
    absolute_pressure_tolerance_m2_s2: f64,
    absolute_flow_tolerance_m2_s: f64,
) -> Result<StokesWindkesselSolution, StokesBoundaryError> {
    let (reduced, equivalent) = reduce_windkessel_ports(
        assembly.system,
        ports,
        loads,
        previous,
        older_history,
        timestep_s,
        absolute_pressure_tolerance_m2_s2,
        absolute_flow_tolerance_m2_s,
    )?;
    let relative_tolerance = assembly.params.tol;
    let prepared = assembly.prepare_pressure_response(ports)?;
    let solved = prepared.solve_resistance_ports(&equivalent, absolute_pressure_tolerance_m2_s2)?;
    let energy = prepared.kinetic_energy(&solved.flow.field.x)?;
    report_windkessel_ports(
        &reduced,
        solved.flow,
        energy,
        solved.pressure_offsets_m2_s2,
        solved.outward_volume_rates_m2_s,
        relative_tolerance,
        absolute_pressure_tolerance_m2_s2,
        absolute_flow_tolerance_m2_s,
    )
}

#[allow(clippy::too_many_arguments)]
fn reduce_windkessel_ports(
    system: &FluxSystem<'_>,
    ports: &[&[usize]],
    loads: &[StokesWindkesselPort],
    previous: StokesWindkesselHistory<'_>,
    older_history: Option<(StokesWindkesselHistory<'_>, f64)>,
    timestep_s: f64,
    absolute_pressure_tolerance_m2_s2: f64,
    absolute_flow_tolerance_m2_s: f64,
) -> Result<(Vec<ReducedWindkesselPort>, Vec<StokesResistancePort>), StokesBoundaryError> {
    let step = StokesBdfStep::new(timestep_s, older_history.map(|(_, dt)| dt)).map_err(invalid)?;
    let valid_history = |history: StokesWindkesselHistory<'_>| {
        history.velocity.len() == system.n_u
            && history.compliance_pressures_m2_s2.len() == ports.len()
            && history
                .velocity
                .iter()
                .chain(history.compliance_pressures_m2_s2)
                .all(|value| value.is_finite())
    };
    if loads.len() != ports.len()
        || !valid_history(previous)
        || older_history.is_some_and(|(history, _)| !valid_history(history))
        || ![
            absolute_pressure_tolerance_m2_s2,
            absolute_flow_tolerance_m2_s,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
    {
        return Err(invalid(
            "complete finite paired fluid/storage histories and tolerances required",
        ));
    }
    if ports.is_empty() || ports.iter().any(|port| port.is_empty()) {
        return Err(invalid("at least one nonempty controlled port required"));
    }
    for &edge in ports.iter().flat_map(|port| port.iter()) {
        if edge >= system.mesh.edges.len() {
            return Err(StokesBoundaryError::EdgeOutOfRange { edge });
        }
        if !system.mesh.is_boundary(edge) {
            return Err(StokesBoundaryError::InteriorEdge { edge });
        }
    }
    let mut reduced = Vec::with_capacity(loads.len());
    let mut equivalent = Vec::with_capacity(loads.len());
    for (i, &load) in loads.iter().enumerate() {
        let previous_q = velocity_port_flow(system, ports[i], previous.velocity);
        let older = older_history.map(|(history, _)| {
            (
                velocity_port_flow(system, ports[i], history.velocity),
                history.compliance_pressures_m2_s2[i],
            )
        });
        let port = ReducedWindkesselPort::new(
            load,
            step,
            previous_q,
            previous.compliance_pressures_m2_s2[i],
            older,
        )
        .map_err(invalid)?;
        equivalent.push(StokesResistancePort {
            reservoir_pressure_m2_s2: port.pressure_offset,
            resistance_s_inv: port.resistance,
        });
        reduced.push(port);
    }
    Ok((reduced, equivalent))
}

#[allow(clippy::too_many_arguments)]
fn report_windkessel_ports<S>(
    reduced: &[ReducedWindkesselPort],
    flow: S,
    kinetic_energy_m4_s2: f64,
    pressure_offsets_m2_s2: Vec<f64>,
    outward_volume_rates_m2_s: Vec<f64>,
    relative_tolerance: f64,
    absolute_pressure_tolerance_m2_s2: f64,
    absolute_flow_tolerance_m2_s: f64,
) -> Result<StokesWindkesselSolution<S>, StokesBoundaryError> {
    if !kinetic_energy_m4_s2.is_finite()
        || kinetic_energy_m4_s2 < 0.0
        || pressure_offsets_m2_s2.len() != reduced.len()
        || outward_volume_rates_m2_s.len() != reduced.len()
    {
        return Err(invalid(
            "unrepresentable hydraulic endpoint energy or incomplete port result",
        ));
    }
    let mut balances = Vec::with_capacity(reduced.len());
    let mut maximum_pressure_residual_m2_s2 = 0.0_f64;
    let mut maximum_flow_residual_m2_s = 0.0_f64;
    for (i, port) in reduced.iter().enumerate() {
        // solve_resistance_ports already re-integrates Q from the accepted BDM
        // field. Reuse that measured value, never a matrix-predicted flow.
        let (balance, pressure_scale, flow_scale) = port
            .report(outward_volume_rates_m2_s[i], pressure_offsets_m2_s2[i])
            .map_err(invalid)?;
        let pressure_tolerance =
            absolute_pressure_tolerance_m2_s2 + relative_tolerance * pressure_scale;
        let flow_tolerance = absolute_flow_tolerance_m2_s + relative_tolerance * flow_scale;
        if !pressure_tolerance.is_finite() || !flow_tolerance.is_finite() {
            return Err(invalid(
                "unrepresentable dynamic hydraulic convergence allowance",
            ));
        }
        let pressure_residual = balance.pressure_law_residual_m2_s2.abs();
        let flow_residual = balance.compliance_flow_residual_m2_s.abs();
        if pressure_residual > pressure_tolerance || flow_residual > flow_tolerance {
            return Err(StokesBoundaryError::WindkesselLawNotConverged {
                port: i,
                pressure_residual_m2_s2: pressure_residual,
                pressure_tolerance_m2_s2: pressure_tolerance,
                flow_residual_m2_s: flow_residual,
                flow_tolerance_m2_s: flow_tolerance,
            });
        }
        maximum_pressure_residual_m2_s2 = maximum_pressure_residual_m2_s2.max(pressure_residual);
        maximum_flow_residual_m2_s = maximum_flow_residual_m2_s.max(flow_residual);
        balances.push(balance);
    }
    Ok(StokesWindkesselSolution {
        flow,
        kinetic_energy_m4_s2,
        pressure_offsets_m2_s2,
        outward_volume_rates_m2_s,
        balances,
        maximum_pressure_residual_m2_s2,
        maximum_flow_residual_m2_s,
    })
}
