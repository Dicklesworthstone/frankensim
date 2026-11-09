//! Fixed-mass heating and free expansion on the same original source isobar.
//!
//! The original Cp owner supplies temperature, enthalpy and sensible entropy;
//! the original engineering-length paths supply free geometry and density.
//! Constant-pressure work is p*Delta V, so Delta U = Delta H - p*Delta V.
//! Cp is never relabeled Cv. This is a uniform quasistatic specimen, with only
//! pressure-volume work, not a constrained solid, pressure solver or moving mesh.

/// Original Cp/conductivity and free density on native implicit thermal networks.
#[path = "expansion/conduction.rs"]
pub mod conduction;

/// Explicit liquid applicability of the same original caloric/expansion path.
#[path = "expansion/liquid.rs"]
pub mod liquid;

use super::IsobaricHeatCapacityLaw;
use crate::thermal::transport::{
    ThermalTransportError,
    expansion::{FreeExpansionDensitySample, PreparedFreeExpansionDensity},
    positive,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::MaterialCard;
use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};

/// Original caloric and free-expansion sources qualified at one maintained pressure.
/// Preparation checks the complete shared material and query context once.
#[derive(Clone, Debug)]
pub struct IsobaricFreeExpansion<'a> {
    caloric: &'a IsobaricHeatCapacityLaw,
    density: &'a PreparedFreeExpansionDensity,
    temperatures_k: Vec<f64>,
    identity: ContentHash,
}

/// One immutable heating/cooling proposal, including its actual free geometry.
/// Energy and entropy are differences; no absolute internal-energy datum is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricExpansionHeatStep {
    /// Original caloric/geometry binding, including the caloric datum.
    pub path_identity: ContentHash,
    /// Original free geometry and density at the starting temperature.
    pub initial: FreeExpansionDensitySample,
    /// Proposed free geometry and density at the accepted heat endpoint.
    pub state: FreeExpansionDensitySample,
    /// Conserved specimen mass [kg].
    pub mass_kg: f64,
    /// Original maintained source pressure [Pa].
    pub pressure_pa: f64,
    /// Signed requested heat into the specimen [J].
    pub supplied_heat_j: f64,
    /// Starting volume from conserved mass and original density [m3].
    pub initial_volume_m3: f64,
    /// Proposed volume from conserved mass and original density [m3].
    pub final_volume_m3: f64,
    /// Integrated original geometry change, retained even when rounded endpoint
    /// volumes or densities coincide. It need not equal their floating-point subtraction.
    pub volume_change_m3: f64,
    /// Direct source integral m*Cp*dT [J].
    pub enthalpy_change_j: f64,
    /// Work delivered to the surroundings, positive for expansion.
    pub boundary_work_j: f64,
    /// Enthalpy change minus pressure-volume work [J].
    pub internal_energy_change_j: f64,
    /// Specimen entropy change from integral m*Cp/T dT, not entropy production.
    pub entropy_change_j_k: f64,
    /// Requested heat minus represented enthalpy change. This is numerical
    /// resolution, not a bound on source uncertainty or equilibration error.
    pub heat_residual_j: f64,
}

impl<'a> IsobaricFreeExpansion<'a> {
    /// Bind original Cp, free lengths and reference density from the SAME card,
    /// phase/process state and exact typed source path. Every selected expansion
    /// and density claim must explicitly support the caloric pressure axis;
    /// an extra pressure coordinate on an otherwise unqualified source is insufficient.
    /// Free expansion and matched-specimen applicability remain caller declarations.
    pub fn try_new(
        card: &MaterialCard,
        caloric: &'a IsobaricHeatCapacityLaw,
        density: &'a PreparedFreeExpansionDensity,
    ) -> Result<Self, ThermalTransportError> {
        let [lower, upper] = caloric.query_endpoints();
        if card.content_hash() != caloric.material_card_identity()
            || card.content_hash() != density.material_card_identity()
            || density.temperature_axis() != caloric.temperature_axis()
            || density.query_endpoints() != [lower, upper]
        {
            return Err(invalid(
                "isobaric free expansion requires one original card and exact source path",
            ));
        }
        for answer in density
            .source_lines()
            .iter()
            .map(|line| &line.knots[0].1)
            .chain(std::iter::once(density.reference_density_answer()))
        {
            let claim = card
                .claims()
                .claim(answer.receipt.selected)
                .ok_or(invalid("original expansion or density claim is missing"))?;
            if !claim
                .validity
                .bounds()
                .contains_key(caloric.pressure_axis())
            {
                return Err(invalid(
                    "free expansion and reference density must support the original isobar",
                ));
            }
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.isobaric-free-expansion.v1");
        hash.update(caloric.identity().as_bytes());
        hash.update(density.identity().as_bytes());
        Ok(Self {
            caloric,
            density,
            temperatures_k: density.temperatures_k().collect(),
            identity: hash.finalize(),
        })
    }

    /// Identity of the retained caloric and free-expansion binding.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Unmodified original Cp law, including its source receipts and datum.
    #[must_use]
    pub const fn caloric(&self) -> &'a IsobaricHeatCapacityLaw {
        self.caloric
    }
    /// Unmodified original free geometry and reference density.
    #[must_use]
    pub const fn density(&self) -> &'a PreparedFreeExpansionDensity {
        self.density
    }

    /// Exact d(1/rho)/dh [m³/J] on this SAME source isobar.
    ///
    /// With rho=rho_ref/J and dh/dT=Cp(T), this is J_T/(rho_ref*Cp).
    /// Cp is evaluated at the original inverse temperature, never frozen at
    /// the density reference or replaced by Cv. Original expansion knots use
    /// their right tangent, except the final endpoint which uses the left.
    /// A zero or negative volume response is retained without clipping.
    pub fn specific_volume_derivative_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<f64, ThermalTransportError> {
        let temperature = self
            .caloric
            .temperature_k(h)
            .map_err(|_| invalid("enthalpy outside the original isobaric expansion path"))?;
        let cp = self.caloric.specific_heat_cp_j_kg_k(temperature)?;
        let tangent = self
            .density
            .expansion()
            .at_temperature_with_tangent(temperature)?;
        let volume_slope = tangent.relative_volume_derivative_per_k();
        self.specific_volume_derivative_from_tangent(volume_slope, cp)
    }

    fn specific_volume_derivative_from_tangent(
        &self,
        volume_slope: f64,
        cp: f64,
    ) -> Result<f64, ThermalTransportError> {
        let reference_density = self.density.reference_density_answer().evidence.value.value;
        let derivative = (volume_slope / reference_density) / cp;
        if !derivative.is_finite() || (volume_slope != 0.0 && derivative == 0.0) {
            return Err(invalid(
                "isobaric specific-volume tangent is unrepresentable",
            ));
        }
        Ok(derivative)
    }

    /// Heat/cool the conserved mass at the source's maintained absolute pressure.
    /// The existing analytic Cp inverse handles every original capacity breakpoint.
    /// Direct interval integrals avoid subtracting large datum-dependent enthalpies.
    /// Original expansion tangents integrate the piecewise-cubic volume exactly
    /// for the admitted affine strain model, up to floating-point roundoff.
    /// A local update searches original breakpoints once and visits only the
    /// crossed intervals; no scan of unrelated parts of the source path is needed.
    /// No source queries, allocations or hashes occur during this evaluation.
    /// Unresolved heat, unsupported temperatures or unrepresentable outputs refuse.
    pub fn advance_heat(
        &self,
        initial_temperature_k: f64,
        mass_kg: f64,
        supplied_heat_j: f64,
    ) -> Result<IsobaricExpansionHeatStep, ThermalTransportError> {
        if !positive(mass_kg) || !supplied_heat_j.is_finite() {
            return Err(invalid(
                "isobaric heating requires positive finite mass and finite heat",
            ));
        }
        let specific_heat = supplied_heat_j / mass_kg;
        if !specific_heat.is_finite() || (supplied_heat_j != 0.0 && specific_heat == 0.0) {
            return Err(invalid(
                "isobaric specific heat increment is unrepresentable",
            ));
        }
        let final_k = self
            .caloric
            .temperature_after_specific_heat(initial_temperature_k, specific_heat)?;
        let dh = self
            .caloric
            .specific_enthalpy_change_j_kg(initial_temperature_k, final_k)?;
        if supplied_heat_j != 0.0
            && (dh == 0.0 || dh.is_sign_positive() != supplied_heat_j.is_sign_positive())
        {
            return Err(invalid(
                "isobaric heat increment is unresolved at the initial temperature",
            ));
        }
        let enthalpy_change_j = product(mass_kg, dh)?;
        let entropy_change_j_k = product(
            mass_kg,
            self.caloric
                .specific_entropy_change_at_temperatures_j_kg_k(initial_temperature_k, final_k)?,
        )?;
        let initial = self.density.at_temperature(initial_temperature_k)?;
        let state = self.density.at_temperature(final_k)?;
        let initial_volume_m3 = mass_kg / initial.density_kg_m3();
        let final_volume_m3 = mass_kg / state.density_kg_m3();
        let reference_volume =
            mass_kg / self.density.reference_density_answer().evidence.value.value;
        if ![initial_volume_m3, final_volume_m3, reference_volume]
            .into_iter()
            .all(positive)
        {
            return Err(invalid("isobaric specimen volume is unrepresentable"));
        }
        let volume_change_m3 = product(
            reference_volume,
            self.relative_volume_change(initial_temperature_k, final_k)?,
        )?;
        let pressure_pa = self.caloric.pressure_pa();
        let boundary_work_j = product(pressure_pa, volume_change_m3)?;
        let internal_energy_change_j = enthalpy_change_j - boundary_work_j;
        let heat_residual_j = supplied_heat_j - enthalpy_change_j;
        if !internal_energy_change_j.is_finite() || !heat_residual_j.is_finite() {
            return Err(invalid("isobaric energy differences are unrepresentable"));
        }
        Ok(IsobaricExpansionHeatStep {
            path_identity: self.identity,
            initial,
            state,
            mass_kg,
            pressure_pa,
            supplied_heat_j,
            initial_volume_m3,
            final_volume_m3,
            volume_change_m3,
            enthalpy_change_j,
            boundary_work_j,
            internal_energy_change_j,
            entropy_change_j_k,
            heat_residual_j,
        })
    }

    fn relative_volume_change(
        &self,
        initial: f64,
        final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        if initial == final_k {
            return Ok(0.0);
        }
        let (low, high, sign) = if initial < final_k {
            (initial, final_k, 1.0)
        } else {
            (final_k, initial, -1.0)
        };
        let geometry = self.density.expansion();
        let mut left = low;
        let mut sum = 0.0;
        let mut correction = 0.0;
        let start = self.temperatures_k.partition_point(|t| *t <= low);
        for right in self.temperatures_k[start..]
            .iter()
            .copied()
            .take_while(|t| *t < high)
            .chain(std::iter::once(high))
        {
            // Every source stretch is affine on this original interval. Expand
            // the cubic product difference, retaining even quadratic changes
            // when opposite principal expansions cancel their linear terms.
            let tangent = geometry.at_temperature_with_tangent(left)?;
            let s = tangent.geometry().principal_stretches();
            let slope = tangent.principal_stretch_derivatives_per_k();
            let dt = right - left;
            let d = [
                product(slope[0], dt)?,
                product(slope[1], dt)?,
                product(slope[2], dt)?,
            ];
            let linear = d[0] * s[1] * s[2] + s[0] * d[1] * s[2] + s[0] * s[1] * d[2];
            let quadratic = d[0] * d[1] * s[2] + d[0] * s[1] * d[2] + s[0] * d[1] * d[2];
            let change = linear + quadratic + d[0] * d[1] * d[2];
            let corrected = change - correction;
            let next = sum + corrected;
            correction = (next - sum) - corrected;
            sum = next;
            left = right;
        }
        if !sum.is_finite() {
            return Err(invalid("original free-volume change is unrepresentable"));
        }
        Ok(sign * sum)
    }
}

impl ThermalLaw for IsobaricFreeExpansion<'_> {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity.as_bytes())
    }

    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::Enthalpy
    }

    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        self.caloric.temperature_k(h)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.caloric.maximum_temperature_slope()
    }

    /// Current material density from the original reference density and J(T).
    /// Stored mass and the reference mesh remain owned by the native timeline.
    /// A geometric consumer explicitly chooses how that volume affects heat
    /// transfer; this does not turn free expansion into constrained-solid data.
    fn bulk_density_kg_m3(&self, h: f64) -> Result<Option<f64>, DepositionError> {
        let temperature = self.caloric.temperature_k(h)?;
        self.density
            .at_temperature(temperature)
            .map(|state| Some(state.density_kg_m3()))
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.caloric
            .specific_entropy_change_j_kg_k(initial_h, final_h)
    }
}

fn invalid(message: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(message)
}
fn product(a: f64, b: f64) -> Result<f64, ThermalTransportError> {
    let value = a * b;
    if !value.is_finite() || (a != 0.0 && b != 0.0 && value == 0.0) {
        Err(invalid(
            "isobaric source increment or extensive product is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}

#[cfg(test)]
#[path = "expansion/tests.rs"]
mod tests;
