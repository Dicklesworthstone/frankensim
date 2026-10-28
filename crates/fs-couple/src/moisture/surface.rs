//! Finite surface exchange in the original dry-basis moisture convention.
//!
//! The source coefficient h_s [m/s] is defined by
//! j_water,in = rho_dry h_s (M_equilibrium - M_surface). It is a surface
//! emission coefficient on dry-basis moisture, not a vapor-concentration or
//! partial-pressure transfer coefficient. The original slab supplies rho_dry;
//! no gas density, equilibrium water inventory or new material evaluator is
//! inferred. Coatings or airflow may motivate source data but are not solved.

use super::*;

/// Source key whose normalization is the dry-basis flux law above.
pub const DRY_BASIS_SURFACE_EMISSION_PROPERTY: &str = "dry_basis_surface_emission_coefficient";
/// Surface emission coefficient dimensions [m/s].
pub const SURFACE_EMISSION_DIMS: Dims = Dims([1, 0, -1, 0, 0, 0]);

/// Pinned, nonnegative surface transfer source at fixed temperature and actual RH.
/// A zero coefficient denotes an impermeable surface. Absence of this object
/// on a sorption law retains the original instantaneous surface equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub struct DryBasisSurfaceTransfer {
    card: MaterialCard,
    source: ResolvedMaterialStatePoint,
}

impl DryBasisSurfaceTransfer {
    /// Admit the named source law. Only typed temperature and relative humidity
    /// are accepted; the coefficient may vary with endpoint RH. Dependence on
    /// evolving surface moisture, hysteresis, vapor pressure or gas velocity is
    /// outside this boundary reduction and cannot be supplied as a hidden axis.
    pub fn new(
        card: MaterialCard,
        point: QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<Self, MoistureError> {
        admit_environment(&point)?;
        let requirement = ScalarPropertyRequirement::try_new(
            DRY_BASIS_SURFACE_EMISSION_PROPERTY,
            SURFACE_EMISSION_DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        let source = resolve_material_state_point(&card, &point, &[requirement], selection)?;
        if !source
            .property(DRY_BASIS_SURFACE_EMISSION_PROPERTY)
            .expect("required surface emission property")
            .has_only_owned_and_fixed_external_conditions(
                "relative_humidity",
                quantity(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            )
        {
            return Err(MoistureError::Invalid(
                "surface exchange needs static source support for the owned RH coordinate",
            ));
        }
        Ok(Self { card, source })
    }

    /// Original immutable material source and selection at the binding point.
    pub const fn source(&self) -> &ResolvedMaterialStatePoint {
        &self.source
    }

    pub(super) fn at(
        &self,
        temperature: f64,
        rh: f64,
    ) -> Result<(f64, ResolvedMaterialStatePoint), MoistureError> {
        if self.source.query().axes()["temperature"].to_bits() != temperature.to_bits()
            || !rh.is_finite()
            || !(0.0..=1.0).contains(&rh)
        {
            return Err(MoistureError::Invalid(
                "surface exchange needs its fixed temperature and actual RH in [0,1]",
            ));
        }
        let point = self
            .source
            .query()
            .with("relative_humidity", rh)
            .map_err(|_| MoistureError::Invalid("invalid relative humidity"))?;
        let source = self.source.at_query_point(&self.card, &point)?;
        let coefficient = source
            .property(DRY_BASIS_SURFACE_EMISSION_PROPERTY)
            .expect("pinned surface emission property")
            .value_si();
        Ok((coefficient, source))
    }
}

pub(super) fn admit_environment(point: &QueryPoint) -> Result<(), MoistureError> {
    if point.axis_quantities().get("temperature")
        != Some(&quantity(QuantityKind::AbsoluteTemperature))
        || point.axis_quantities().get("relative_humidity")
            != Some(&quantity(QuantityKind::Moisture(
                MoistureBasis::RelativeHumidity,
            )))
        || !point
            .axes()
            .get("temperature")
            .copied()
            .is_some_and(positive)
        || !point
            .axes()
            .get("relative_humidity")
            .is_some_and(|v| v.is_finite() && (0.0..=1.0).contains(v))
        || point
            .axes()
            .keys()
            .any(|key| key != "temperature" && key != "relative_humidity")
    {
        return Err(MoistureError::Invalid(
            "surface laws require only typed temperature and relative humidity coordinates",
        ));
    }
    Ok(())
}

#[derive(Default)]
pub(super) struct ResolvedMoistureFace {
    pub equilibrium: Option<f64>,
    pub coefficient_m_s: Option<f64>,
    pub sorption: Option<ResolvedMaterialStatePoint>,
    pub transfer: Option<ResolvedMaterialStatePoint>,
}

impl ResolvedMoistureFace {
    /// Series combination of half-cell diffusion and finite surface resistance.
    /// Returns dt*k_eff/dx and the fraction of (M_eq-M_cell) at the surface.
    pub(super) fn coupling(
        &self,
        dt: f64,
        dx: f64,
        diffusivity: f64,
        scale: f64,
    ) -> Result<(f64, f64), MoistureError> {
        if self.equilibrium.is_none() {
            return Ok((0.0, 0.0));
        }
        let Some(h) = self.coefficient_m_s else {
            // Preserve the original Dirichlet path and its arithmetic exactly.
            return Ok((2.0 * scale * diffusivity, 1.0));
        };
        if h == 0.0 {
            return Ok((0.0, 0.0));
        }
        let diffusion_speed = 2.0 * diffusivity / dx;
        if !positive(diffusion_speed) || !positive(h) {
            return Err(MoistureError::Invalid(
                "surface and half-cell transfer speeds must be representable",
            ));
        }
        let lo = h.min(diffusion_speed);
        let hi = h.max(diffusion_speed);
        let effective = lo / (1.0 + lo / hi);
        let coefficient = dt / dx * effective;
        let fraction = if h <= diffusion_speed {
            let ratio = h / diffusion_speed;
            ratio / (1.0 + ratio)
        } else {
            1.0 / (1.0 + diffusion_speed / h)
        };
        if !positive(coefficient) || !positive(fraction) {
            return Err(MoistureError::Invalid(
                "finite surface transfer coefficient overflow or underflow",
            ));
        }
        Ok((coefficient, fraction))
    }
}
