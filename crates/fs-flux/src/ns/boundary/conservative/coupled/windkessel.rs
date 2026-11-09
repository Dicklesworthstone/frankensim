//! Transported-density BDF1 momentum with the original finite RCR port storage.

use super::{
    CoupledDensityMomentumStep, DensityMomentumControls, DensityTransportCell, FluxSystem,
    MaterialFluxNumerics, MaterialNavierStokesBoundary, StokesRheologyControls, invalid, poll,
};
use crate::material::mixed::{physical_windkessel, with_kinematic_windkessel};
use crate::material::{
    MaterialFluxError, MaterialStokesWindkesselHistory, MaterialStokesWindkesselPort,
    MaterialStokesWindkesselSolution,
};
use crate::ns::boundary::StokesBoundaryError;
use crate::ns::boundary::conservative::DensityMomentumError;
use crate::ns::boundary::flow_control::{
    FlowControl, PortLaw, reduce_windkessel_ports, report_windkessel_ports,
};

/// Original refusal from either coupled material transport or the hydraulic load.
#[derive(Debug, Clone, PartialEq)]
pub enum DensityMomentumWindkesselError {
    /// Density/momentum admission, nonlinear solve, observation or cancellation.
    Momentum(DensityMomentumError),
    /// Physical port parameters, normalization, original RCR laws or reporting.
    Hydraulic(MaterialFluxError),
}

impl core::fmt::Display for DensityMomentumWindkesselError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Momentum(error) => error.fmt(f),
            Self::Hydraulic(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for DensityMomentumWindkesselError {}
impl From<DensityMomentumError> for DensityMomentumWindkesselError {
    fn from(error: DensityMomentumError) -> Self {
        Self::Momentum(error)
    }
}
impl From<MaterialFluxError> for DensityMomentumWindkesselError {
    fn from(error: MaterialFluxError) -> Self {
        Self::Hydraulic(error)
    }
}
impl From<StokesBoundaryError> for DensityMomentumWindkesselError {
    fn from(error: StokesBoundaryError) -> Self {
        Self::Hydraulic(MaterialFluxError::Boundary(error))
    }
}

type Result<T> = core::result::Result<T, DensityMomentumWindkesselError>;

impl FluxSystem<'_> {
    /// Advance transported density, BDF1 momentum and finite RCR storage together.
    ///
    /// The original passive load obeys `p=p_res+Rp*Q+L*D_t(Q)+x` and
    /// `C*D_t(x)=Q-x/Rd`, with outward volume flow Q and compliance pressure x
    /// relative to its reservoir. Both use the SAME BDF1 interval as momentum.
    /// The native owner reduces these laws once from immutable old velocity and
    /// compliance. Each density/viscosity Picard trial then resolves the load's
    /// pressure through its own native pressure responses. Rejected iterations
    /// never advance the external storage or transport another interval.
    ///
    /// The accepted endpoint is reassembled at its actual transported density
    /// and velocity. Original RCR laws, storage, resistor losses and BDF defect
    /// are evaluated using its measured port flow and added pressure. The same
    /// pressure increments remain in the returned field and boundary reactions.
    /// All traction ports may carry loads: their reservoir laws supply pressure
    /// references without an extra uncontrolled datum.
    ///
    /// Existing SI port/history types apply: R and L are per-depth coefficients,
    /// C is volume-per-depth compliance, and histories are physical Pa. Parameters,
    /// port order and compliance references stay fixed through the interval;
    /// reservoir pressure increments are held over it. Flow is m²/s, and the
    /// returned external energies/powers remain J/m and W/m. Multiply these by
    /// controls.depth_m for whole-domain J/W, matching the density step's already
    /// depth-integrated inventories and kinetic energy. Reference density only
    /// conditions the equations; actual cell density supplies fluid inertia.
    ///
    /// This advances a lumped mechanical compliance, not an external material,
    /// gas EOS, thermal inventory or moving reservoir geometry. Resistor loss
    /// and numerical BDF defect remain separate from transported caloric data.
    /// Incoming density/inventory and momentum donors remain caller-owned under
    /// Open/OpenReservoir semantics, including reversal at any loaded port.
    ///
    /// All reconstruction, CFL, divergence and storage limits of the original
    /// endpoint-advector step apply. Callbacks must be pure. Cancellation is
    /// checked around native solves/reporting and throughout density transport;
    /// no fluid or external-storage endpoint escapes any complete-step refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn density_momentum_coupled_bdf1_windkessel_step<const N: usize>(
        &self,
        numerics: MaterialFluxNumerics,
        dynamic_viscosity: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> f64,
        force_density: &dyn Fn(usize, [f64; 2], &DensityTransportCell<N>) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        ports: &[&[usize]],
        loads: &[MaterialStokesWindkesselPort],
        previous: MaterialStokesWindkesselHistory<'_>,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        controls: DensityMomentumControls,
        nonlinear: StokesRheologyControls,
        absolute_pressure_tolerance_pa: f64,
        absolute_flow_tolerance_m2_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialStokesWindkesselSolution<CoupledDensityMomentumStep<N>>> {
        poll(&mut cancelled)?;
        let reference = controls.reference_density_kg_m3;
        // The nested result preserves both original error owners without
        // replacing transport/cancellation failures by a normalization error.
        let native = with_kinematic_windkessel(
            reference,
            ports,
            loads,
            previous,
            None,
            absolute_pressure_tolerance_pa,
            absolute_flow_tolerance_m2_s,
            |loads, history, _| {
                Ok((|| -> Result<_> {
                    poll(&mut cancelled)?;
                    let pressure_tolerance = absolute_pressure_tolerance_pa / reference;
                    let (reduced, equivalent) = reduce_windkessel_ports(
                        self,
                        ports,
                        loads,
                        history,
                        None,
                        timestep_s,
                        pressure_tolerance,
                        absolute_flow_tolerance_m2_s,
                    )?;
                    poll(&mut cancelled)?;
                    let solved = self.coupled_density_momentum_bdf1(
                        numerics,
                        dynamic_viscosity,
                        force_density,
                        boundaries,
                        history.velocity,
                        cells,
                        inflow,
                        timestep_s,
                        controls,
                        nonlinear,
                        Some(FlowControl {
                            ports,
                            law: PortLaw::Resistance {
                                loads: &equivalent,
                                absolute_pressure_tolerance_m2_s2: pressure_tolerance,
                            },
                        }),
                        &mut cancelled,
                    )?;
                    // The uniform reference is NEVER a replacement density
                    // metric. Recover native per-depth units from the accepted
                    // heterogeneous-density observation, not kinetic_energy(u).
                    let physical_energy = solved.flow.step.observations.kinetic_energy_after_j;
                    let energy_per_depth = physical_energy / controls.depth_m;
                    let native_energy = energy_per_depth / reference;
                    if !energy_per_depth.is_finite()
                        || !native_energy.is_finite()
                        || (physical_energy > 0.0 && native_energy == 0.0)
                    {
                        return Err(
                            invalid("unrepresentable physical hydraulic fluid energy").into()
                        );
                    }
                    poll(&mut cancelled)?;
                    let reported = report_windkessel_ports(
                        &reduced,
                        solved.flow,
                        native_energy,
                        solved.pressure_offsets_m2_s2,
                        solved.outward_volume_rates_m2_s,
                        numerics.linear_tolerance,
                        pressure_tolerance,
                        absolute_flow_tolerance_m2_s,
                    )?;
                    poll(&mut cancelled)?;
                    Ok(reported)
                })())
            },
        )??;
        let physical = physical_windkessel(native, reference, |flow, _| Ok(flow))?;
        poll(&mut cancelled)?;
        Ok(physical)
    }
}
