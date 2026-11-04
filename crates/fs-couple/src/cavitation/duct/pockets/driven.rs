//! Prescribed endpoint inflow through the original DG pressure/work port.
//!
//! The held flow drives the SAME liquid, contact events and gas-energy states.
//! It is a source boundary, not a pump curve, valve law or moving-mesh model.
//! Positive flow is into the duct at either endpoint. Piecewise-held histories
//! may change flow between accepted intervals; no acoustic state is reset.

use super::*;
use caloric::CaloricGasObservation;

/// Accepted source exchange, separate from every internal cavity/gas exchange.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctDriveStep {
    pub port: PortId,
    pub inflow_m3_s: f64,
    /// Signed prescribed inflow times the actual interval duration [m3].
    pub transported_volume_m3: f64,
    pub end_absolute_pressure_pa: f64,
    /// Integral of the original pressure PERTURBATION at this port [Pa s].
    pub acoustic_pressure_impulse_pa_s: f64,
    /// Integral p' Q dt into original acoustic storage [J]. This, together
    /// with gas heat, is the external input to the coupled energy balance.
    pub acoustic_work_j: f64,
    /// p_reference * integral Q dt [J]. This background-pressure transport is
    /// reported separately, not added to the small-strain acoustic Hamiltonian
    /// or deposited as heat. Adding it to acoustic_work_j gives absolute p Q.
    pub reference_pressure_transport_work_j: f64,
}

pub(super) struct Drive {
    pub port: PortId,
    pub port_index: usize,
    area_m2: f64,
    /// B_pocket^T b_drive, in the existing pocket order.
    pub pressure_rates: [f64; MAX_DUCT_POCKETS],
}

impl Drive {
    pub(super) fn new<F: CavitationProperties + FluidAcoustics + Clone>(
        duct: &DgDuct<DgFluidMaterial<F>>,
        port: PortId,
        sites: &[Site],
    ) -> Result<Self, DuctCavityError> {
        let area_m2 = if duct.inlet() == DgBoundary::Port(port) {
            duct.inlet_area_m2()
        } else if duct.outlet() == DgBoundary::Port(port) {
            duct.outlet_area_m2()
        } else {
            return Err(DuctCavityError::InvalidInput(
                "flow drive must name an original endpoint port",
            ));
        };
        let view = duct.continuous_operator();
        let port_index = view
            .ports()
            .iter()
            .position(|(id, _)| *id == port)
            .ok_or(DuctCavityError::InvalidInput("undeclared flow-drive port"))?;
        let mut pressure_rates = [0.0; MAX_DUCT_POCKETS];
        for (i, site) in sites.iter().enumerate() {
            pressure_rates[i] = dot(
                &view.ports()[site.port_index].1,
                &view.ports()[port_index].1,
            );
        }
        if !positive(area_m2) || pressure_rates.iter().any(|r| !r.is_finite()) {
            return Err(DuctCavityError::InvalidInput(
                "flow-drive geometry or pressure rates overflow",
            ));
        }
        Ok(Self {
            port,
            port_index,
            area_m2,
            pressure_rates,
        })
    }

    pub(super) fn admit(
        &self,
        inflow: f64,
        dt: f64,
        sound_speed: f64,
        limits: DuctCavityLimits,
    ) -> Result<(), DuctCavityError> {
        let volume = inflow * dt;
        let reference_work = limits.reference_pressure_pa * volume;
        if !inflow.is_finite()
            || (inflow / self.area_m2 / sound_speed).abs() > limits.maximum_interface_mach
            || !volume.is_finite()
            || (inflow != 0.0 && volume == 0.0)
            || !reference_work.is_finite()
            || self
                .pressure_rates
                .iter()
                .any(|r| !(r * inflow).is_finite())
        {
            return Err(DuctCavityError::OutsideDomain(
                "prescribed inflow exceeds native Mach or represented volume/work bounds",
            ));
        }
        Ok(())
    }

    pub(super) fn report(
        &self,
        inflow: f64,
        dt: f64,
        reference: f64,
        endpoint_pressure: f64,
        impulse: f64,
        work: f64,
    ) -> Result<DuctDriveStep, DuctCavityError> {
        let volume = inflow * dt;
        let absolute_pressure = reference + endpoint_pressure;
        let reference_work = reference * volume;
        if ![volume, absolute_pressure, impulse, work, reference_work]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(DuctCavityError::OutsideDomain(
                "accepted flow-drive exchange is unrepresentable",
            ));
        }
        Ok(DuctDriveStep {
            port: self.port,
            inflow_m3_s: inflow,
            transported_volume_m3: volume,
            end_absolute_pressure_pa: absolute_pressure,
            acoustic_pressure_impulse_pa_s: impulse,
            acoustic_work_j: work,
            reference_pressure_transport_work_j: reference_work,
        })
    }
}

impl<F: CavitationProperties + FluidAcoustics + Clone> CavitatingDuctPockets<F> {
    /// Declare one original endpoint as a prescribed-flow source. Every other
    /// original port must have exactly one pocket. The drive cannot also be a
    /// pocket; all original material, wall memory and local volume bounds stay
    /// active. Construction admits the initial state at zero source flow.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_driven_port(
        duct: DgDuct<DgFluidMaterial<F>>,
        initial: DgState,
        time_s: f64,
        specifications: &[DuctPocket],
        gas_pockets: &[TrappedGasPocket],
        caloric_pockets: &[CaloricGasPocket],
        driven_port: PortId,
        limits: DuctCavityLimits,
    ) -> Result<Self, DuctCavityError> {
        Self::bind(
            duct,
            initial,
            time_s,
            specifications,
            gas_pockets,
            caloric_pockets,
            Some(driven_port),
            limits,
        )
    }

    /// Advance with one signed inflow held over this interval and zero gas heat.
    /// Positive inflow compresses the liquid column or displaces its pockets;
    /// both endpoints use the native INTO-duct convention. A change from the
    /// previous interval is included in initial contact resolution. The source
    /// cannot silently create a cavity at its own undeclared liquid endpoint.
    pub fn advance_with_port_flow<Cancel: FnMut() -> bool>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        held_inflow_m3_s: f64,
        cancelled: &mut Cancel,
    ) -> Result<DuctPocketsStep, DuctCavityError> {
        self.advance_with_port_flow_and_gas_heat(
            end_time_s,
            numerics,
            held_inflow_m3_s,
            &|_, _| Ok(0.0),
            cancelled,
        )
    }

    /// Mechanical and thermal input on the same accepted native trajectory.
    /// The heat callback has the same deterministic/rejected-trial semantics
    /// as `advance_with_gas_heat`. All gas/liquid state and receipts publish
    /// together after source, event, work and cancellation admission.
    pub fn advance_with_port_flow_and_gas_heat<Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        held_inflow_m3_s: f64,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<DuctPocketsStep, DuctCavityError>
    where
        Heat: Fn(f64, &CaloricGasObservation) -> Result<f64, DuctCavityError>,
        Cancel: FnMut() -> bool,
    {
        if self.drive.is_none() {
            return Err(DuctCavityError::InvalidInput(
                "no endpoint flow drive was declared",
            ));
        }
        self.advance_kernel(end_time_s, numerics, heat, held_inflow_m3_s, cancelled)
    }

    pub(super) fn drive_pressure_rate(&self, pocket: usize, inflow: f64) -> f64 {
        self.drive
            .as_ref()
            .map_or(0.0, |d| d.pressure_rates[pocket] * inflow)
    }
}
