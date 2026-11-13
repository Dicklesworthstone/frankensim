//! Hydraulic endpoint feedback through the original DG pressure/work port.
//!
//! Held flow and a smooth Q(time, absolute pressure) law drive the SAME liquid,
//! contact events and gas-energy states. Every stage and contact root uses its
//! actual pressure, and transported volume and p'Q work are integrated on that
//! same trajectory. Positive flow is into the duct at either endpoint.
//! The passive reservoir law is a declared quasi-steady resistance, not a pump
//! curve or an inertial valve model. No boundary pressure is clipped or reset.

use super::*;
use caloric::CaloricGasObservation;

/// Accepted source exchange, separate from every internal cavity/gas exchange.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctDriveStep {
    pub port: PortId,
    /// Instantaneous endpoint inflow [m3/s]. Constant for a held-flow request.
    pub inflow_m3_s: f64,
    pub initial_inflow_m3_s: f64,
    /// Actual signed volume integral divided by the physical duration [m3/s].
    pub average_inflow_m3_s: f64,
    /// Integral of actual signed inflow on the native trajectory [m3].
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
    /// Present only when this SAME advance used the passive reservoir law.
    pub reservoir: Option<DuctReservoirBalance>,
}

/// Constant-pressure reservoir of the same frozen liquid, through an authored
/// positive hydraulic resistance. A calibrated resistance or a separately
/// justified quasi-steady constitutive reduction must supply R [Pa s/m3].
/// No geometry, pump head, loss coefficient or material datum is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctReservoirResistance {
    absolute_pressure_pa: f64,
    resistance_pa_s_m3: f64,
}

impl DuctReservoirResistance {
    pub fn try_new(
        absolute_pressure_pa: f64,
        resistance_pa_s_m3: f64,
    ) -> Result<Self, DuctCavityError> {
        if !positive(absolute_pressure_pa) || !positive(resistance_pa_s_m3) {
            return Err(DuctCavityError::InvalidInput(
                "positive finite reservoir pressure and hydraulic resistance required",
            ));
        }
        Ok(Self {
            absolute_pressure_pa,
            resistance_pa_s_m3,
        })
    }

    pub const fn absolute_pressure_pa(self) -> f64 {
        self.absolute_pressure_pa
    }
    pub const fn resistance_pa_s_m3(self) -> f64 {
        self.resistance_pa_s_m3
    }

    /// Signed Q=(p_res-p_port)/R. Finite tensile port pressures are not clipped;
    /// the original liquid/cavity owner independently admits that field state.
    pub fn inflow_m3_s(self, absolute_port_pressure_pa: f64) -> Result<f64, DuctCavityError> {
        let drop = self.absolute_pressure_pa - absolute_port_pressure_pa;
        let flow = drop / self.resistance_pa_s_m3;
        if !absolute_port_pressure_pa.is_finite()
            || !drop.is_finite()
            || !flow.is_finite()
            || (drop != 0.0 && flow == 0.0)
        {
            return Err(DuctCavityError::OutsideDomain(
                "reservoir pressure/flow is unrepresentable",
            ));
        }
        Ok(flow)
    }
}

/// Separate energy ports of the actual passive reservoir-driven trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctReservoirBalance {
    pub absolute_pressure_pa: f64,
    pub resistance_pa_s_m3: f64,
    /// p_res integral Q dt, signed work supplied by the external reservoir [J].
    pub supplied_absolute_work_j: f64,
    /// (p_res-p_reference) integral Q dt in the liquid's energy chart [J].
    pub supplied_reduced_work_j: f64,
    /// Actual integral R Q^2 dt [J] on the SAME native trajectory. This is
    /// external resistance loss; it is not silently deposited in the liquid.
    pub dissipated_energy_j: f64,
    /// Supplied reduced work minus acoustic work and resistance loss [J].
    pub energy_balance_residual_j: f64,
}

pub(super) type FlowLaw<'a> = dyn Fn(f64, f64) -> Result<f64, DuctCavityError> + 'a;

#[derive(Clone, Copy)]
pub(super) enum DriveSignal<'a> {
    Held(f64),
    Feedback(&'a FlowLaw<'a>),
    Reservoir(DuctReservoirResistance),
}

impl DriveSignal<'_> {
    pub(super) fn sample(self, time_s: f64, pressure_pa: f64) -> Result<f64, DuctCavityError> {
        match self {
            Self::Held(flow) => Ok(flow),
            Self::Feedback(law) => law(time_s, pressure_pa),
            Self::Reservoir(reservoir) => reservoir.inflow_m3_s(pressure_pa),
        }
    }

    pub(super) fn balance(
        self,
        reference_pa: f64,
        volume_m3: f64,
        acoustic_work_j: f64,
        dissipated_energy_j: f64,
        numerics: &DuctCavityNumerics,
    ) -> Result<Option<DuctReservoirBalance>, DuctCavityError> {
        let Self::Reservoir(reservoir) = self else {
            return Ok(None);
        };
        let absolute = reservoir.absolute_pressure_pa * volume_m3;
        let reduced = (reservoir.absolute_pressure_pa - reference_pa) * volume_m3;
        // Cancel the common absolute pressure before arithmetic. Resistance
        // loss has its own integral, so subtraction cannot erase a small loss.
        let residual = reduced - acoustic_work_j - dissipated_energy_j;
        let tolerance = numerics.absolute_energy_tolerance_j
            + numerics.relative_tolerance
                * (reduced.abs() + acoustic_work_j.abs() + dissipated_energy_j.abs());
        if ![absolute, reduced, dissipated_energy_j, residual, tolerance]
            .into_iter()
            .all(f64::is_finite)
            || dissipated_energy_j < 0.0
            || residual.abs() > tolerance
        {
            return Err(DuctCavityError::EnergyBalance);
        }
        Ok(Some(DuctReservoirBalance {
            absolute_pressure_pa: reservoir.absolute_pressure_pa,
            resistance_pa_s_m3: reservoir.resistance_pa_s_m3,
            supplied_absolute_work_j: absolute,
            supplied_reduced_work_j: reduced,
            dissipated_energy_j,
            energy_balance_residual_j: residual,
        }))
    }
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
                "hydraulic inflow exceeds native Mach or represented volume/work bounds",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn report(
        &self,
        initial_inflow: f64,
        final_inflow: f64,
        integrated_volume: f64,
        signal: DriveSignal<'_>,
        dt: f64,
        reference: f64,
        endpoint_pressure: f64,
        impulse: f64,
        work: f64,
        dissipated_energy: f64,
        numerics: &DuctCavityNumerics,
    ) -> Result<DuctDriveStep, DuctCavityError> {
        // The held-flow volume integral is analytically exact. A general law
        // must retain its integrated history, not endpoint flow times dt.
        let volume = match signal {
            DriveSignal::Held(inflow) => inflow * dt,
            _ => integrated_volume,
        };
        let average = volume / dt;
        let absolute_pressure = reference + endpoint_pressure;
        let reference_work = reference * volume;
        if ![
            volume,
            average,
            absolute_pressure,
            impulse,
            work,
            reference_work,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(DuctCavityError::OutsideDomain(
                "accepted flow-drive exchange is unrepresentable",
            ));
        }
        Ok(DuctDriveStep {
            port: self.port,
            inflow_m3_s: final_inflow,
            initial_inflow_m3_s: initial_inflow,
            average_inflow_m3_s: average,
            transported_volume_m3: volume,
            end_absolute_pressure_pa: absolute_pressure,
            acoustic_pressure_impulse_pa_s: impulse,
            acoustic_work_j: work,
            reference_pressure_transport_work_j: reference_work,
            reservoir: signal.balance(reference, volume, work, dissipated_energy, numerics)?,
        })
    }
}

impl<F: CavitationProperties + FluidAcoustics + Clone> CavitatingDuctPockets<F> {
    /// Declare one original endpoint as a hydraulic source. Every other
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
        self.advance_kernel(
            end_time_s,
            numerics,
            heat,
            DriveSignal::Held(held_inflow_m3_s),
            cancelled,
        )
    }

    /// Evolve a smooth, deterministic Q(t, p_absolute) boundary on the SAME
    /// dense liquid/cavity trajectory. The callback is evaluated at every
    /// physical stage and contact root, never just once per requested window.
    /// It must be bounded, side-effect-free and valid at rejected trial states.
    /// Resolve its bandwidth with maximum_step_s and split known discontinuities.
    /// Every current flow is independently checked against the native Mach cap.
    pub fn advance_with_port_flow_law<Flow, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        flow: &Flow,
        cancelled: &mut Cancel,
    ) -> Result<DuctPocketsStep, DuctCavityError>
    where
        Flow: Fn(f64, f64) -> Result<f64, DuctCavityError>,
        Cancel: FnMut() -> bool,
    {
        self.advance_with_port_flow_law_and_gas_heat(
            end_time_s,
            numerics,
            flow,
            &|_, _| Ok(0.0),
            cancelled,
        )
    }

    /// Pressure-dependent hydraulic flow and finite gas heat share one accepted
    /// trajectory, numerical budget and atomic publication boundary.
    pub fn advance_with_port_flow_law_and_gas_heat<Flow, Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        flow: &Flow,
        heat: &Heat,
        cancelled: &mut Cancel,
    ) -> Result<DuctPocketsStep, DuctCavityError>
    where
        Flow: Fn(f64, f64) -> Result<f64, DuctCavityError>,
        Heat: Fn(f64, &CaloricGasObservation) -> Result<f64, DuctCavityError>,
        Cancel: FnMut() -> bool,
    {
        if self.drive.is_none() {
            return Err(DuctCavityError::InvalidInput(
                "no endpoint flow drive was declared",
            ));
        }
        self.advance_kernel(
            end_time_s,
            numerics,
            heat,
            DriveSignal::Feedback(flow),
            cancelled,
        )
    }

    /// A constant-pressure reservoir of this liquid feeds or drains the duct
    /// through its explicitly supplied positive quasi-steady resistance.
    /// Its pressure must remain at or above saturation and within the same frozen
    /// liquid's declared small-strain domain.
    /// Reverse flow is retained. The actual reservoir work and dissipated
    /// energy are admitted before any native liquid/gas state is published.
    pub fn advance_with_reservoir<Cancel: FnMut() -> bool>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        reservoir: DuctReservoirResistance,
        cancelled: &mut Cancel,
    ) -> Result<DuctPocketsStep, DuctCavityError> {
        self.advance_with_reservoir_and_gas_heat(
            end_time_s,
            numerics,
            reservoir,
            &|_, _| Ok(0.0),
            cancelled,
        )
    }

    pub fn advance_with_reservoir_and_gas_heat<Heat, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &DuctCavityNumerics,
        reservoir: DuctReservoirResistance,
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
        if reservoir.absolute_pressure_pa < self.material.saturation_pressure_pa() {
            return Err(DuctCavityError::OutsideDomain(
                "the external liquid reservoir must be at or above its vapor pressure",
            ));
        }
        let reservoir_strain = (reservoir.absolute_pressure_pa - self.limits.reference_pressure_pa)
            / self.duct.acoustic_material().adiabatic_bulk_modulus_pa();
        if !reservoir_strain.is_finite()
            || reservoir_strain.abs() > self.limits.maximum_liquid_strain
        {
            return Err(DuctCavityError::OutsideDomain(
                "reservoir pressure exceeds the frozen liquid's strain domain",
            ));
        }
        self.advance_kernel(
            end_time_s,
            numerics,
            heat,
            DriveSignal::Reservoir(reservoir),
            cancelled,
        )
    }

    pub(super) fn sample_drive(
        &self,
        signal: DriveSignal<'_>,
        time_s: f64,
        coordinates: &[f64],
        coordinate_scale: f64,
        duration_s: f64,
    ) -> Result<f64, DuctCavityError> {
        let Some(drive) = &self.drive else {
            return Ok(0.0);
        };
        // Preserve the held-flow lane's original work and adaptive state size;
        // its analytically known flow does not need another pressure dot product.
        if let DriveSignal::Held(inflow) = signal {
            drive.admit(
                inflow,
                duration_s,
                self.material.liquid_sound_speed_m_s(),
                self.limits,
            )?;
            return Ok(inflow);
        }
        let view = self.duct.continuous_operator();
        let pressure = self.limits.reference_pressure_pa
            + coordinate_scale * dot(&view.ports()[drive.port_index].1, coordinates);
        if !time_s.is_finite() || !pressure.is_finite() {
            return Err(DuctCavityError::OutsideDomain(
                "hydraulic stage pressure/time is nonfinite",
            ));
        }
        let inflow = signal.sample(time_s, pressure)?;
        drive.admit(
            inflow,
            duration_s,
            self.material.liquid_sound_speed_m_s(),
            self.limits,
        )?;
        Ok(inflow)
    }

    pub(super) fn drive_pressure_rate(&self, pocket: usize, inflow: f64) -> f64 {
        self.drive
            .as_ref()
            .map_or(0.0, |d| d.pressure_rates[pocket] * inflow)
    }
}
