//! Original conductivity integrals across a growing deposit and its contact.
//!
//! This adapter binds mass-owned coating geometry to the existing
//! [`ReservoirLayerGeometry`] material/film owner. It integrates the original
//! positive `k(T)` chart across the ACTUAL coating temperature drop and solves
//! the contact-face temperature; conductivity is not frozen at a bulk trace.
//! The separate-temperature deposition kernel retains both caloric inventories,
//! molecular exchange, time integration and the equal/opposite interface heat.
//!
//! The layer is uniform, isotropic, source-free and quasi-steady through its
//! thickness. Its mass and constant-density morphology are frozen during one
//! constitutive interval; rebind this allocation-free object as mass evolves.
//! No source fit, pore model, extra storage or resolved temperature field is
//! inferred. A source-qualified owner may lend its `conductivity_curve()`;
//! borrowing a curve alone does not establish material/caloric identity.

use crate::thermal::conduction::{
    ConductivityCurve, LayerConductionConfig, ReservoirLayerGeometry, ReservoirLayerResponse,
};
use fs_thermal::enthalpy::DepositionError;
use fs_thermal::phase_geometry::FrontGeometry;
use fs_thermal::sublimation::SublimationError;
use fs_thermal::sublimation::coating::resistance::{
    ResistiveInterfaceLaw, geometry::ConductiveCoatingGeometry,
};

/// A borrowed original conductivity chart and one actual positive coating mass.
///
/// Geometry preparation and successful flux evaluations allocate nothing and
/// perform no material-card query. The original layer solver owns at most
/// `config.max_iterations` supported scalar trials per heat-rate query; each
/// original conductivity operation visits at most its bounded knot count.
#[derive(Clone, Copy)]
pub struct ConductiveCoatingInterface<'a> {
    geometry: ConductiveCoatingGeometry,
    front: FrontGeometry,
    conductivity: &'a ConductivityCurve,
    layer: ReservoirLayerGeometry,
    contact_resistance_k_w: f64,
    config: LayerConductionConfig,
    checkpoint: Option<&'a dyn Fn() -> Result<(), DepositionError>>,
}

impl<'a> ConductiveCoatingInterface<'a> {
    /// Bind actual coating mass, the original k(T) curve and inner contact area.
    ///
    /// The unit-conductivity query extracts only the ORIGINAL geometric metric:
    /// its conductance is `C = 1/integral(ds/A)` in metres. No constant material
    /// conductivity enters the flux calculation. The original geometry retains
    /// positive sub-ULP shell thicknesses without subtracting rounded radii.
    /// The contact's areal resistance is divided by the fixed INNER area once.
    /// Bare/extinct films and invalid/unrepresentable metrics refuse.
    pub fn new(
        geometry: ConductiveCoatingGeometry,
        deposit_mass_kg: f64,
        conductivity: &'a ConductivityCurve,
        contact_resistance_m2_k_w: f64,
        config: LayerConductionConfig,
    ) -> Result<Self, DepositionError> {
        if !contact_resistance_m2_k_w.is_finite()
            || contact_resistance_m2_k_w < 0.0
            || config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(DepositionError::InvalidInput(
                "integrated coating needs finite nonnegative contact resistance and bounded flux controls",
            ));
        }
        let metric = geometry
            .at_mass(deposit_mass_kg, 1.0, 0.0)
            .map_err(geometry_error)?;
        let layer = ReservoirLayerGeometry::from_geometric_multiplier(
            metric.conductance_w_k,
            metric.contact_area_m2,
        )?;
        let contact_resistance_k_w = contact_resistance_m2_k_w / metric.contact_area_m2;
        if !contact_resistance_k_w.is_finite()
            || (contact_resistance_m2_k_w > 0.0 && contact_resistance_k_w == 0.0)
        {
            return Err(DepositionError::Arithmetic);
        }
        Ok(Self {
            geometry,
            front: metric.front,
            conductivity,
            layer,
            contact_resistance_k_w,
            config,
            checkpoint: None,
        })
    }

    /// Retain an immutable scope cancellation gate inside every original layer
    /// iteration, including calls through [`ResistiveInterfaceLaw`]. Without
    /// this gate a trait call remains bounded by the explicit layer budget and
    /// the outer deposition owner polls around each constitutive evaluation.
    #[must_use]
    pub fn with_checkpoint<'scope>(
        self,
        checkpoint: &'scope dyn Fn() -> Result<(), DepositionError>,
    ) -> ConductiveCoatingInterface<'scope>
    where
        'a: 'scope,
    {
        ConductiveCoatingInterface {
            geometry: self.geometry,
            front: self.front,
            conductivity: self.conductivity,
            layer: self.layer,
            contact_resistance_k_w: self.contact_resistance_k_w,
            config: self.config,
            checkpoint: Some(checkpoint),
        }
    }

    #[must_use]
    pub const fn geometry(&self) -> ConductiveCoatingGeometry {
        self.geometry
    }

    /// Deposit-only mass, thickness and outer area from the original geometry.
    #[must_use]
    pub const fn front(&self) -> FrontGeometry {
        self.front
    }

    /// The borrowed original knots, not samples fitted at trial temperatures.
    #[must_use]
    pub const fn conductivity_curve(&self) -> &'a ConductivityCurve {
        self.conductivity
    }

    #[must_use]
    pub const fn contact_resistance_k_w(&self) -> f64 {
        self.contact_resistance_k_w
    }

    /// Original nonlinear layer/contact response in [support, deposit] order.
    ///
    /// Positive power leaves the substrate and enters the deposit. The reported
    /// material-face temperature is the coating's inner contact trace. A finite
    /// contact jump may keep that trace inside source support even when the
    /// substrate temperature is outside it; the substrate is never sampled as
    /// coating material. Perfect contact requires both traces to be supported.
    /// Original typed source/iteration/arithmetic/cancellation errors survive.
    /// The layer's heat-continuity allowance is separate from the outer thermal
    /// energy and adaptive time tolerances; it is not a physical source-fit bound.
    pub fn evaluate(
        &self,
        temperatures_k: [f64; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ReservoirLayerResponse, DepositionError> {
        let mut check = || {
            if let Some(scope) = self.checkpoint {
                scope()?;
            }
            checkpoint()
        };
        // The original material-to-bath evaluator keeps source support on the
        // coating side of the contact. Reverse its ordered transfer exactly
        // once to express support -> deposit, including both derivatives.
        let mut response = self.layer.evaluate_material_to_bath(
            self.conductivity,
            [temperatures_k[1], temperatures_k[0]],
            self.contact_resistance_k_w,
            self.config,
            &mut check,
        )?;
        response.power_to_reservoir_w = -response.power_to_reservoir_w;
        response.heat_residual_w = -response.heat_residual_w;
        let [deposit, support] = response.endpoint_derivative_w_k;
        response.endpoint_derivative_w_k = [-support, -deposit];
        Ok(response)
    }
}

impl ResistiveInterfaceLaw for ConductiveCoatingInterface<'_> {
    fn heat_rate_w(
        &self,
        substrate_temperature_k: f64,
        deposit_temperature_k: f64,
    ) -> Result<f64, SublimationError> {
        self.evaluate([substrate_temperature_k, deposit_temperature_k], || Ok(()))
            .map(|response| response.power_to_reservoir_w)
            .map_err(interface_error)
    }
}

fn geometry_error(error: SublimationError) -> DepositionError {
    match error {
        SublimationError::Invalid(message) => DepositionError::InvalidInput(message),
        SublimationError::TemperatureOutOfRange { .. } => DepositionError::MaterialDomain,
        SublimationError::NonConverged => DepositionError::Budget("coating geometry"),
        SublimationError::Cancelled => DepositionError::Cancelled,
    }
}

// The scalar thermal trait intentionally has a small low-level error vocabulary.
// `evaluate` above remains available when the caller needs the original detail.
fn interface_error(error: DepositionError) -> SublimationError {
    match error {
        DepositionError::Cancelled => SublimationError::Cancelled,
        DepositionError::Budget(_) | DepositionError::EnergyBalance(_) => {
            SublimationError::NonConverged
        }
        DepositionError::InvalidInput(message) => SublimationError::Invalid(message),
        DepositionError::MaterialDomain => {
            SublimationError::Invalid("coating interface left original conductivity source support")
        }
        DepositionError::Arithmetic => {
            SublimationError::Invalid("source-integrated coating heat is unrepresentable")
        }
    }
}

#[cfg(test)]
#[path = "conductive/tests.rs"]
mod tests;
