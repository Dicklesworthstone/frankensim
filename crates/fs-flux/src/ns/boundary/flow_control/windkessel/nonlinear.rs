//! Dynamic external loads inside the original convection/constitutive iteration.

use super::super::{FlowControl, PortLaw};
use super::{
    StokesBoundaryError, StokesWindkesselHistory, StokesWindkesselPort, StokesWindkesselSolution,
    reduce_windkessel_ports, report_windkessel_ports,
};
use crate::ns::boundary::{NavierStokesBoundary, NavierStokesControls, NavierStokesStep};
use crate::ns::{FluxParams, FluxSystem};

impl FluxSystem<'_> {
    /// Advance laminar convection, source rheology and external RCR/inertance
    /// at the SAME fully implicit BDF endpoint.
    ///
    /// Every original Picard iteration updates convection and the supplied
    /// `law(cell, point, shear_rate)` [m²/s], then solves all external pressure
    /// laws through that current native mobility. The endpoint is reassembled
    /// at its accepted shear and velocity; both original nonlinear equations
    /// and integrated port laws must pass. Pressure increments relax together
    /// with velocity and fluid pressure, preserving the same coupled load.
    ///
    /// BDF histories, compliance coordinates and units follow the linear
    /// Windkessel entry. Native natural/dissipative backflow and incoming
    /// reservoir momentum remain explicit in `boundaries`. The RCR report
    /// accounts only for its external storage/losses: convected kinetic energy,
    /// viscous heat and backflow-policy work are distinct fluid terms. No
    /// arbitrary nonlinear passivity or turbulent closure is implied.
    ///
    /// # Errors
    /// Preserves source, history, boundary and original nonlinear-budget
    /// refusals, plus both unreduced dynamic-law checks. No partial fluid or
    /// capacitor state is returned and both histories remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[StokesWindkesselPort],
        previous: StokesWindkesselHistory<'_>,
        older_history: Option<(StokesWindkesselHistory<'_>, f64)>,
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_pressure_tolerance_m2_s2: f64,
        absolute_flow_tolerance_m2_s: f64,
    ) -> Result<StokesWindkesselSolution<NavierStokesStep>, StokesBoundaryError> {
        self.bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports_checked(
            params,
            law,
            force,
            boundaries,
            ports,
            loads,
            previous,
            older_history,
            timestep_s,
            controls,
            absolute_pressure_tolerance_m2_s2,
            absolute_flow_tolerance_m2_s,
            &|| true,
        )
    }

    /// The same fully implicit convective/storage step with constant viscosity.
    /// `older_history=None` uses BDF1 startup; no artificial convection freeze
    /// or prescribed boundary profile is introduced.
    ///
    /// # Errors
    /// Retains all nonlinear dynamic-port checks and requires positive finite nu.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_variable_step_windkessel_ports(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[StokesWindkesselPort],
        previous: StokesWindkesselHistory<'_>,
        older_history: Option<(StokesWindkesselHistory<'_>, f64)>,
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_pressure_tolerance_m2_s2: f64,
        absolute_flow_tolerance_m2_s: f64,
    ) -> Result<StokesWindkesselSolution<NavierStokesStep>, StokesBoundaryError> {
        if !params.nu.is_finite() || params.nu <= 0.0 {
            return Err(StokesBoundaryError::InvalidInput(
                "invalid Navier–Stokes viscosity",
            ));
        }
        self.bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports(
            params,
            &|_, _, _| params.nu,
            force,
            boundaries,
            ports,
            loads,
            previous,
            older_history,
            timestep_s,
            controls,
            absolute_pressure_tolerance_m2_s2,
            absolute_flow_tolerance_m2_s,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn bdf2_navier_stokes_mixed_rheology_variable_step_windkessel_ports_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[StokesWindkesselPort],
        previous: StokesWindkesselHistory<'_>,
        older_history: Option<(StokesWindkesselHistory<'_>, f64)>,
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_pressure_tolerance_m2_s2: f64,
        absolute_flow_tolerance_m2_s: f64,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<StokesWindkesselSolution<NavierStokesStep>, StokesBoundaryError> {
        let (reduced, equivalent) = reduce_windkessel_ports(
            self,
            ports,
            loads,
            previous,
            older_history,
            timestep_s,
            absolute_pressure_tolerance_m2_s2,
            absolute_flow_tolerance_m2_s,
        )?;
        let solved = self.bdf2_navier_stokes_mixed_prepared_rheology_with_flow_control(
            params,
            &|velocity| {
                let rates = self.shear_rates(velocity);
                if rates.iter().any(|rate| !rate.is_finite()) {
                    return Err(StokesBoundaryError::InvalidInput("nonfinite shear rate"));
                }
                Ok(rates)
            },
            law,
            force,
            boundaries,
            previous.velocity,
            older_history.map(|(history, dt)| (history.velocity, dt)),
            timestep_s,
            controls,
            source_is_valid,
            Some(FlowControl {
                ports,
                law: PortLaw::Resistance {
                    loads: &equivalent,
                    absolute_pressure_tolerance_m2_s2,
                },
            }),
        )?;
        let energy = self.kinetic_energy(&solved.flow.solution.field.x);
        report_windkessel_ports(
            &reduced,
            solved.flow,
            energy,
            solved.pressure_offsets_m2_s2,
            solved.outward_volume_rates_m2_s,
            params.tol,
            absolute_pressure_tolerance_m2_s2,
            absolute_flow_tolerance_m2_s,
        )
    }
}
