//! Constant mass-flow cooling on the original native conducting surface.
//!
//! The bulk stream owns its velocity: U = mdot / (rho_bulk A_flow). Film
//! density still supplies the selected correlation's Reynolds number; it must
//! not replace rho_bulk in the flow conversion. No fluid store, pressure solve,
//! fan law, second heat integrator or fitted material is introduced.

use super::evaluator::ExteriorFluid;
use super::*;

/// Current bulk velocity and the original coefficient/surface tangent.
#[derive(Debug)]
pub struct MassFlowConvectionResponse {
    pub mass_flow_kg_s: f64,
    /// Flow cross section, distinct from the native boundary's heated area.
    pub flow_area_m2: f64,
    /// Original source at the physical bath temperature, not the film state.
    pub bulk_transport: FluidTransportState,
    pub speed_m_s: f64,
    pub convection: ExteriorConvectionResponse,
    /// Absolute dh/dmdot at fixed actual surface/bulk temperatures, flow area,
    /// geometry, pressure and composition. The input is a nonnegative flow
    /// MAGNITUDE, not a signed branch rate or a total fan-system derivative.
    /// `None` means the original correlation supplies no local two-sided Re
    /// partial, including its validity endpoints; it never means zero. Section
    /// laws use their section coefficient's partial, not a downstream mean.
    /// This is a local formula derivative, not a material/model error bound.
    pub coefficient_mass_flow_derivative_w_m2_k_s_kg: Option<f64>,
}

/// Prescribed mass flow through an original gas or liquid exterior law.
///
/// Borrow a prepared `ImplicitGasConvection` or `ImplicitLiquidConvection`;
/// all original source selections, caloric identity, temperature support,
/// property conventions, surface-viscosity corrections and physical geometry
/// remain attached. The wrapped law's authored speed is replaced by the bulk
/// velocity at EACH actual bath temperature. Only forced laws are admitted:
/// a scalar flow magnitude cannot specify natural/mixed flow orientation.
///
/// `mass_flow_kg_s` is nonnegative and `flow_area_m2` is the stream's positive
/// cross section. It is not the thermal patch area. Direction and the physical
/// leading/heated entrance remain those of the retained law; this adapter does
/// not infer a new entrance after flow reversal. A zero rate is never raised to
/// an artificial speed floor and still needs the original correlation's domain.
///
/// During native surface elimination the bath, mass rate and area are fixed,
/// hence dU/dTs = 0. The existing film-property or wall-viscosity chain supplies
/// the exact same surface tangent as a correctly rebound speed law. Across bath
/// updates the original bulk density changes the speed. No bath-temperature,
/// hydraulic-network or complete coupled-system Jacobian is claimed.
pub struct PrescribedMassFlowConvection<'a> {
    source: &'a dyn ExteriorFluid,
    mass_flow_kg_s: f64,
    flow_area_m2: f64,
    mass_flux_kg_m2_s: f64,
}

impl<'a> PrescribedMassFlowConvection<'a> {
    pub fn for_gas(
        source: &'a ImplicitGasConvection<'_, '_>,
        mass_flow_kg_s: f64,
        flow_area_m2: f64,
    ) -> Result<Self, ExteriorConvectionError> {
        Self::new(source, mass_flow_kg_s, flow_area_m2)
    }

    pub fn for_liquid(
        source: &'a ImplicitLiquidConvection<'_, '_>,
        mass_flow_kg_s: f64,
        flow_area_m2: f64,
    ) -> Result<Self, ExteriorConvectionError> {
        Self::new(source, mass_flow_kg_s, flow_area_m2)
    }

    fn new(
        source: &'a dyn ExteriorFluid,
        mass_flow_kg_s: f64,
        flow_area_m2: f64,
    ) -> Result<Self, ExteriorConvectionError> {
        if !mass_flow_kg_s.is_finite()
            || mass_flow_kg_s < 0.0
            || !flow_area_m2.is_finite()
            || flow_area_m2 <= 0.0
        {
            return Err(ExteriorConvectionError::Invalid(
                "prescribed convection needs a finite nonnegative mass-flow magnitude and positive flow area",
            ));
        }
        // Admit geometry class without sampling the material or pretending
        // that an arbitrary reference temperature is the actual bulk state.
        with_speed(source.settings().law, 0.0)?;
        let mass_flux_kg_m2_s = mass_flow_kg_s / flow_area_m2;
        if !mass_flux_kg_m2_s.is_finite() || (mass_flow_kg_s > 0.0 && mass_flux_kg_m2_s <= 0.0) {
            return Err(ExteriorConvectionError::Invalid(
                "prescribed mass flux is unrepresentable",
            ));
        }
        Ok(Self {
            source,
            mass_flow_kg_s,
            flow_area_m2,
            mass_flux_kg_m2_s,
        })
    }

    pub const fn mass_flow_kg_s(&self) -> f64 {
        self.mass_flow_kg_s
    }

    pub const fn flow_area_m2(&self) -> f64 {
        self.flow_area_m2
    }

    /// Resolve the original bulk and correlation states without mutating them.
    /// A local rate partial is optional; its absence does not invalidate an
    /// otherwise admitted held-rate native surface solve.
    pub fn evaluate(
        &self,
        surface_temperature_k: f64,
        bath_temperature_k: f64,
    ) -> Result<MassFlowConvectionResponse, ExteriorConvectionError> {
        let bound = self.bind(bath_temperature_k)?;
        let convection = evaluator::evaluate(&bound, surface_temperature_k, bath_temperature_k)?;
        let coefficient_mass_flow_derivative_w_m2_k_s_kg =
            self.mass_flow_partial(&bound, &convection)?;
        Ok(MassFlowConvectionResponse {
            mass_flow_kg_s: self.mass_flow_kg_s,
            flow_area_m2: self.flow_area_m2,
            bulk_transport: bound.bulk_transport,
            speed_m_s: bound.speed_m_s,
            convection,
            coefficient_mass_flow_derivative_w_m2_k_s_kg,
        })
    }

    fn bind(&self, bath: f64) -> Result<BoundFlow<'_>, ExteriorConvectionError> {
        let bulk_transport = self.source.sample_transport(bath, false)?.state;
        let speed_m_s = self.mass_flux_kg_m2_s / bulk_transport.density_kg_m3();
        if !speed_m_s.is_finite() || (self.mass_flow_kg_s > 0.0 && speed_m_s <= 0.0) {
            return Err(ExteriorConvectionError::Invalid(
                "actual bulk-density convection velocity is unrepresentable",
            ));
        }
        let mut settings = self.source.settings();
        settings.law = with_speed(settings.law, speed_m_s)?;
        Ok(BoundFlow {
            source: self.source,
            settings,
            bulk_temperature_k: bath,
            bulk_transport,
            speed_m_s,
        })
    }

    fn mass_flow_partial(
        &self,
        bound: &BoundFlow<'_>,
        response: &ExteriorConvectionResponse,
    ) -> Result<Option<f64>, ExteriorConvectionError> {
        // At zero rate the magnitude map has no admitted two-sided derivative.
        // The primal still goes through the unchanged original correlation.
        if self.mass_flow_kg_s == 0.0 {
            return Ok(None);
        }
        let length = bound.settings.characteristic_length;
        let conductivity =
            ThermalConductivity::new(response.transport.thermal_conductivity_w_m_k());
        let reynolds = response.nusselt.groups()["Re"];
        let prandtl = response.transport.prandtl();
        let dh_dre = match bound.settings.law {
            SourceCellConvectionLaw::ForcedPlateSection { .. } => {
                response.coefficient.sensitivity.d_qoi.get("Re").copied()
            }
            SourceCellConvectionLaw::ForcedTubeSection { start, end, .. } => {
                evaluate_turbulent_circular_duct_section(
                    reynolds,
                    prandtl,
                    start.value() / length.value(),
                    end.value() / length.value(),
                )?
                .heat_transfer_coefficient_with_sensitivities(conductivity, length)?
                .sensitivity
                .d_qoi
                .get("Re")
                .copied()
            }
            SourceCellConvectionLaw::ForcedLaminarTubeSection { start, end, .. } => {
                evaluate_circular_duct_section(
                    reynolds,
                    prandtl,
                    start.value() / length.value(),
                    end.value() / length.value(),
                )?
                .heat_transfer_coefficient_with_sensitivities(conductivity, length)?
                .sensitivity
                .d_qoi
                .get("Re")
                .copied()
            }
            SourceCellConvectionLaw::Forced { .. }
            | SourceCellConvectionLaw::ForcedSphere { .. } => response
                .nusselt
                .evidence()
                .sensitivity
                .d_qoi
                .get("Re")
                .map(|partial| partial * conductivity.value() / length.value()),
            _ => unreachable!("constructor admits only supported forced geometry"),
        };
        dh_dre
            .map(|partial| {
                // A source-owned constant-Nu direction remains exactly zero
                // even if an independently formed Re/mdot would overflow.
                if partial == 0.0 {
                    Ok(0.0)
                } else {
                    finite(partial * (reynolds / self.mass_flow_kg_s))
                }
            })
            .transpose()
    }
}

fn with_speed(
    law: SourceCellConvectionLaw,
    speed_m_s: f64,
) -> Result<SourceCellConvectionLaw, ExteriorConvectionError> {
    Ok(match law {
        SourceCellConvectionLaw::Forced {
            correlation,
            settings,
            ..
        } => SourceCellConvectionLaw::Forced {
            correlation,
            speed_m_s,
            settings,
        },
        SourceCellConvectionLaw::ForcedPlateSection {
            correlation, start, ..
        } => SourceCellConvectionLaw::ForcedPlateSection {
            correlation,
            speed_m_s,
            start,
        },
        SourceCellConvectionLaw::ForcedTubeSection { start, end, .. } => {
            SourceCellConvectionLaw::ForcedTubeSection {
                speed_m_s,
                start,
                end,
            }
        }
        SourceCellConvectionLaw::ForcedLaminarTubeSection { start, end, .. } => {
            SourceCellConvectionLaw::ForcedLaminarTubeSection {
                speed_m_s,
                start,
                end,
            }
        }
        SourceCellConvectionLaw::ForcedSphere { .. } => {
            SourceCellConvectionLaw::ForcedSphere { speed_m_s }
        }
        _ => {
            return Err(ExteriorConvectionError::Invalid(
                "prescribed mass-flow exterior requires a forced pointwise, plate-strip or circular-tube law",
            ));
        }
    })
}

// Only this ephemeral view changes velocity. All material/tangent queries and
// native caloric identity retain the original prepared owner. Bulk-held laws
// reuse the already resolved bulk snapshot instead of querying it twice.
struct BoundFlow<'a> {
    source: &'a dyn ExteriorFluid,
    settings: ExteriorConvectionSettings,
    bulk_temperature_k: f64,
    bulk_transport: FluidTransportState,
    speed_m_s: f64,
}

impl ExteriorFluid for BoundFlow<'_> {
    fn caloric(&self) -> &dyn ThermalLaw {
        self.source.caloric()
    }
    fn energy_bounds(&self) -> [f64; 2] {
        self.source.energy_bounds()
    }
    fn settings(&self) -> ExteriorConvectionSettings {
        self.settings
    }
    fn sample_transport(
        &self,
        temperature: f64,
        derivatives: bool,
    ) -> Result<TransportSample, ExteriorConvectionError> {
        if !derivatives && temperature == self.bulk_temperature_k {
            Ok(TransportSample {
                state: self.bulk_transport,
                derivatives: [0.0; 4],
            })
        } else {
            self.source.sample_transport(temperature, derivatives)
        }
    }
    fn wall_viscosity(
        &self,
        temperature: f64,
    ) -> Result<(f64, f64, bool), ExteriorConvectionError> {
        self.source.wall_viscosity(temperature)
    }
    fn expansion(
        &self,
        temperature: f64,
        sample: &TransportSample,
        derivative: bool,
    ) -> Result<(f64, f64), ExteriorConvectionError> {
        self.source.expansion(temperature, sample, derivative)
    }
}

impl ImplicitGrayExteriorLaw for PrescribedMassFlowConvection<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.source.caloric()
    }
    fn energy_bounds_j_kg(&self) -> [f64; 2] {
        self.source.energy_bounds()
    }
    fn surface_temperature_bounds_k(&self) -> [f64; 2] {
        self.source.settings().surface_temperature_bounds_k
    }
    fn flux_state_at_surface(
        &self,
        energy: f64,
        surface: f64,
        bath: f64,
    ) -> Result<ImplicitGrayExteriorFluxState, DepositionError> {
        let bound = self
            .bind(bath)
            .map_err(|_| DepositionError::MaterialDomain)?;
        evaluator::flux_state_at_surface(&bound, energy, surface, bath)
    }
    fn state_at_surface(
        &self,
        energy: f64,
        surface: f64,
        bath: f64,
    ) -> Result<ImplicitGrayExteriorState, DepositionError> {
        let bound = self
            .bind(bath)
            .map_err(|_| DepositionError::MaterialDomain)?;
        evaluator::state_at_surface(&bound, energy, surface, bath)
    }
}
