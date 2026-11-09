//! Bulk/wall material transport on the original bare tube-bank correlation.
/// Finite liquid throughflow with source-exact calorics and refreshed bank properties.
pub mod stream;

use super::*;
use fs_convection::tube_bank::{self, TubeBankError, TubeBankEvaluation, TubeBankGeometry};

#[derive(Debug)]
pub enum MaterialTubeBankError {
    Material(MaterialError),
    Correlation(TubeBankError),
}
impl core::fmt::Display for MaterialTubeBankError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "{error}"),
            Self::Correlation(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for MaterialTubeBankError {}
impl From<MaterialError> for MaterialTubeBankError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<TubeBankError> for MaterialTubeBankError {
    fn from(error: TubeBankError) -> Self {
        Self::Correlation(error)
    }
}

/// Both source owners and exact snapshots used by the bank mean and wall correction.
#[derive(Debug)]
pub struct MaterialTubeBankConvection<'a, B: FluidTransport + ?Sized, S: FluidTransport + ?Sized> {
    bulk: &'a B,
    surface: &'a S,
    bulk_transport: FluidTransportState,
    surface_transport: FluidTransportState,
    evaluation: TubeBankEvaluation,
    coefficient: Evidence<HeatTransferCoefficient>,
}
impl<B: FluidTransport + ?Sized, S: FluidTransport + ?Sized> MaterialTubeBankConvection<'_, B, S> {
    pub const fn bulk_fluid(&self) -> &B {
        self.bulk
    }
    pub const fn surface_fluid(&self) -> &S {
        self.surface
    }
    pub const fn bulk_transport(&self) -> &FluidTransportState {
        &self.bulk_transport
    }
    pub const fn surface_transport(&self) -> &FluidTransportState {
        &self.surface_transport
    }
    pub const fn evaluation(&self) -> &TubeBankEvaluation {
        &self.evaluation
    }
    pub const fn coefficient(&self) -> &Evidence<HeatTransferCoefficient> {
        &self.coefficient
    }
    /// Same bulk conductivity and outside diameter; no second material query.
    /// The reference is an effective bulk temperature appropriate for the
    /// bank mean, not automatically the inlet temperature of a warming stream.
    pub fn robin_boundary(
        &self,
        reference: Temperature,
    ) -> Result<CorrelationRobinBoundary, TubeBankError> {
        self.evaluation.robin_boundary(
            ThermalConductivity::new(self.bulk_transport.thermal_conductivity_w_m_k()),
            reference,
        )
    }
}

/// Mean bank convection using UPSTREAM unobstructed approach speed [m/s].
/// Geometry supplies the maximum inter-tube speed used by the correlation.
/// Bulk rho/mu define Re; bulk mu/Cp/k define Pr and h; wall-fluid mu/Cp/k
/// supply Pr_surface. Each material is sampled exactly once and retained.
///
/// Both owners must describe the same fluid/composition/pressure at their
/// declared mean-bulk and wall temperatures. This generic transport seam cannot
/// establish that identity. No air preset, wall-solid conductivity, fan model,
/// finite-fluid energy balance or short-pin/fin correction is inferred.
pub fn evaluate_tube_bank_from_material<
    'a,
    B: FluidTransport + ?Sized,
    S: FluidTransport + ?Sized,
>(
    bulk: &'a B,
    surface: &'a S,
    approach_speed_m_s: f64,
    geometry: TubeBankGeometry,
) -> Result<MaterialTubeBankConvection<'a, B, S>, MaterialTubeBankError> {
    let bulk_transport = bulk.transport_state()?;
    let surface_transport = surface.transport_state()?;
    let reynolds = bulk_transport.reynolds(approach_speed_m_s, geometry.diameter().value())?;
    let evaluation = tube_bank::evaluate(
        geometry,
        reynolds,
        bulk_transport.prandtl(),
        surface_transport.prandtl(),
    )?;
    let coefficient = evaluation.heat_transfer_coefficient(ThermalConductivity::new(
        bulk_transport.thermal_conductivity_w_m_k(),
    ))?;
    Ok(MaterialTubeBankConvection {
        bulk,
        surface,
        bulk_transport,
        surface_transport,
        evaluation,
        coefficient,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_convection::tube_bank::TubeBankLayout;

    #[test]
    fn g1_material_bank_keeps_bulk_and_wall_prandtl_and_lowers_actual_coefficient() {
        let bulk = FluidTransportState::try_new(1000.0, 3500.0, 0.5, 0.001).unwrap();
        let surface = FluidTransportState::try_new(970.0, 4000.0, 0.8, 0.0007).unwrap();
        let geometry = TubeBankGeometry::try_new(
            TubeBankLayout::Inline,
            4,
            Length::new(0.01),
            Length::new(0.02),
            Length::new(0.02),
        )
        .unwrap();
        let result = evaluate_tube_bank_from_material(&bulk, &surface, 1.0, geometry).unwrap();
        assert!((result.evaluation().groups()["Re_max"] - 20_000.0).abs() < 1e-10);
        assert!((result.evaluation().groups()["Pr"] - 7.0).abs() < 1e-14);
        assert!((result.evaluation().groups()["Pr_surface"] - 3.5).abs() < 1e-14);
        assert!((result.coefficient().value.value() - 50.0 * 298.36506027896674).abs() < 1e-8);
        assert!(core::ptr::eq(result.bulk_fluid(), &bulk));
        assert!(core::ptr::eq(result.surface_fluid(), &surface));
        assert_eq!(*result.bulk_transport(), bulk);
        assert_eq!(*result.surface_transport(), surface);
        assert_eq!(
            result
                .robin_boundary(Temperature::new(300.0))
                .unwrap()
                .coefficient(),
            result.coefficient()
        );
        assert!(evaluate_tube_bank_from_material(&bulk, &surface, 0.0, geometry).is_err());
    }
}
