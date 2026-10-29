//! Native volumetric caloric data on the existing conservative thermal timeline.
//!
//! The timeline stores `m*h`. This adapter uses a caller-declared, fixed reference
//! mass/volume ratio rho0 and `h(T) = integral C_vol(T) dT / rho0`, so a compatible
//! cell stores exactly `V0 * integral C_vol(T) dT`. rho0 is a scenario input, not a
//! newly measured density or a temperature-dependent density law. The caller
//! selects the source's caloric constraint and justifies fixed volume and mass.
//! The existing timeline owns all state, time integration, sources and balances.
//! No pressure work, phase change, moisture exchange or moving volume is added.

use super::{PreparedVolumetricThermalTransport, ThermalTransportError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionSimulation, ThermalLaw, ThermalLink,
    series_contact_conductance_w_k,
};

/// Borrowed source model expressed in the timeline's fixed reference mass basis.
///
/// Zero specific storage is at the explicit datum temperature. Density does not
/// change with temperature; it only relates the caller's mass and volume. Every
/// cell using this law must satisfy that same relation. Use `admit_cell` before
/// inserting cells; sourced links repeat that admission at each evaluation.
#[derive(Clone, Debug)]
pub struct FixedVolumeThermalLaw<'p> {
    path: &'p PreparedVolumetricThermalTransport,
    reference_density: f64,
    datum_temperature: f64,
    specific_storage_span: [f64; 2],
    maximum_slope: f64,
    sensible_entropy: bool,
    identity: ContentHash,
}

impl<'p> FixedVolumeThermalLaw<'p> {
    pub fn try_new(
        path: &'p PreparedVolumetricThermalTransport,
        reference_density_kg_m3: f64,
        datum_temperature_k: f64,
    ) -> Result<Self, ThermalTransportError> {
        if !positive(reference_density_kg_m3) {
            return Err(invalid(
                "fixed reference mass per volume must be finite and positive",
            ));
        }
        path.at_temperature(datum_temperature_k)?;
        let span = path.temperature_span_k();
        let specific_storage_span = [
            path.storage_change_j_m3(datum_temperature_k, span[0])? / reference_density_kg_m3,
            path.storage_change_j_m3(datum_temperature_k, span[1])? / reference_density_kg_m3,
        ];
        if !specific_storage_span.iter().all(|h| h.is_finite())
            || specific_storage_span[0] >= specific_storage_span[1]
            || (span[0] != datum_temperature_k && specific_storage_span[0] == 0.0)
            || (span[1] != datum_temperature_k && specific_storage_span[1] == 0.0)
        {
            return Err(invalid(
                "specific storage is not representable in this reference mass basis",
            ));
        }
        // A positive affine capacity reaches its minimum at a retained knot.
        // Round the ratio upward so the explicit owner's CFL bound is not
        // underestimated by this division's final rounding.
        let mut minimum_capacity = f64::INFINITY;
        for temperature in path.temperatures_k() {
            minimum_capacity = minimum_capacity.min(
                path.at_temperature(temperature)?
                    .volumetric_heat_capacity_j_m3_k(),
            );
        }
        let maximum_slope = (reference_density_kg_m3 / minimum_capacity).next_up();
        if !positive(maximum_slope) {
            return Err(invalid(
                "temperature slope is not representable in this reference mass basis",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.fixed-volume-thermal-law.v1");
        hash.update(path.identity().as_bytes());
        hash.update(&reference_density_kg_m3.to_bits().to_le_bytes());
        hash.update(&datum_temperature_k.to_bits().to_le_bytes());
        Ok(Self {
            path,
            reference_density: reference_density_kg_m3,
            datum_temperature: datum_temperature_k,
            specific_storage_span,
            maximum_slope,
            sensible_entropy: false,
            identity: hash.finalize(),
        })
    }

    /// Explicitly select reversible sensible thermal entropy on this same curve.
    ///
    /// The caller declares the capacity appropriate to the fixed caloric/work
    /// constraint, composition and reference volume. This enables the existing
    /// timeline's relative entropy account; it does not infer an absolute entropy
    /// datum, source entropy, missing work coordinates or an experimental entropy
    /// measurement. Default laws continue to report entropy as unavailable.
    /// The selected interpretation changes the binding identity exactly once.
    #[must_use]
    pub fn with_sensible_entropy(mut self) -> Self {
        if !self.sensible_entropy {
            let mut hash =
                DomainHasher::new("org.frankensim.fs-material.fixed-volume-sensible-entropy.v1");
            hash.update(self.identity.as_bytes());
            self.identity = hash.finalize();
            self.sensible_entropy = true;
        }
        self
    }

    #[must_use]
    pub const fn sensible_entropy_enabled(&self) -> bool {
        self.sensible_entropy
    }

    /// The unchanged source path, including both original lines and receipts.
    #[must_use]
    pub const fn path(&self) -> &'p PreparedVolumetricThermalTransport {
        self.path
    }

    /// Identity includes the source path, mass basis, datum and entropy opt-in.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    #[must_use]
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density
    }

    #[must_use]
    pub const fn datum_temperature_k(&self) -> f64 {
        self.datum_temperature
    }

    #[must_use]
    pub const fn specific_storage_span_j_kg(&self) -> [f64; 2] {
        self.specific_storage_span
    }

    /// Derived specific storage for this fixed mass basis; not a source Cp claim.
    pub fn specific_storage_j_kg(&self, temperature_k: f64) -> Result<f64, ThermalTransportError> {
        let storage = self
            .path
            .storage_change_j_m3(self.datum_temperature, temperature_k)?
            / self.reference_density;
        if !storage.is_finite() || (temperature_k != self.datum_temperature && storage == 0.0) {
            return Err(invalid("specific storage is not representable"));
        }
        Ok(storage)
    }

    /// Delivered mass compatible with a declared fixed reference volume.
    pub fn reference_mass_kg(&self, volume_m3: f64) -> Result<f64, ThermalTransportError> {
        let mass = self.reference_density * volume_m3;
        if !positive(volume_m3) || !positive(mass) {
            return Err(invalid(
                "reference volume and resulting mass must be finite and positive",
            ));
        }
        Ok(mass)
    }

    /// Admit the actual cell's mass/volume basis and supported stored coordinate.
    /// The 64-epsilon comparison allows arithmetic roundoff, not physical density
    /// uncertainty or a material tolerance. This does not certify spatial resolution.
    pub fn admit_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        let expected = self
            .reference_mass_kg(cell.geometry.volume_m3)
            .map_err(|_| DepositionError::MaterialDomain)?;
        if !positive(cell.mass_kg)
            || (cell.mass_kg - expected).abs() > 64.0 * f64::EPSILON * cell.mass_kg.max(expected)
        {
            return Err(DepositionError::InvalidInput(
                "cell mass and reference volume do not match its volumetric thermal law",
            ));
        }
        self.temperature_k(cell.specific_h_j_kg)?;
        Ok(())
    }
}

impl ThermalLaw for FixedVolumeThermalLaw<'_> {
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        let [low, high] = self.specific_storage_span;
        if !specific_h_j_kg.is_finite() || specific_h_j_kg < low || specific_h_j_kg > high {
            return Err(DepositionError::MaterialDomain);
        }
        let span = self.path.temperature_span_k();
        // Preserve the exact original domain endpoints across the mass conversion.
        if specific_h_j_kg == low {
            return Ok(span[0]);
        }
        if specific_h_j_kg == high {
            return Ok(span[1]);
        }
        let storage = specific_h_j_kg * self.reference_density;
        if !storage.is_finite() || (specific_h_j_kg != 0.0 && storage == 0.0) {
            return Err(DepositionError::MaterialDomain);
        }
        self.path
            .temperature_after_storage_change(self.datum_temperature, storage)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        if !self.sensible_entropy {
            return Ok(None);
        }
        let initial_k = self.temperature_k(initial_h_j_kg)?;
        let final_k = self.temperature_k(final_h_j_kg)?;
        if initial_h_j_kg == final_h_j_kg {
            return Ok(Some(0.0));
        }
        let change_h = final_h_j_kg - initial_h_j_kg;
        let entropy = if initial_k == final_k {
            // Distinct stored h coordinates can round to the same T. Retain the
            // actual small dh using the isothermal limit instead of erasing it.
            change_h / initial_k
        } else {
            let entropy_v = self
                .path
                .sensible_entropy_change_j_m3_k(initial_k, final_k)
                .map_err(|_| DepositionError::MaterialDomain)?;
            let storage_v = self
                .path
                .storage_change_j_m3(initial_k, final_k)
                .map_err(|_| DepositionError::MaterialDomain)?;
            // ds/dh is the capacity-weighted mean of 1/T. Use the actual
            // stored dh so inversion roundoff cannot erase a small heat input.
            // In exact arithmetic this is entropy_v/reference_density.
            change_h * (entropy_v / storage_v)
        };
        if !entropy.is_finite()
            || entropy == 0.0
            || entropy.is_sign_positive() != (final_h_j_kg > initial_h_j_kg)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(entropy))
    }
}

/// One existing thermal link evaluated from its two actual material states.
///
/// The caller declares the path's scalar k to be aligned with each transfer leg.
/// Half-cell distances, area and contact resistance are explicit geometry/model
/// inputs. Each leg uses k at its cell's uniform temperature, a finite-volume
/// approximation rather than an exact integral through a nonuniform half-cell.
/// Source preparation is amortized; evaluations query no card and hash nothing.
#[derive(Clone, Debug)]
pub struct FixedVolumeConductionLink<'l, 'p> {
    id: u64,
    cells: [CellId; 2],
    laws: [&'l FixedVolumeThermalLaw<'p>; 2],
    distances_m: [f64; 2],
    area_m2: f64,
    contact_resistance_m2_k_w: f64,
}

impl<'l, 'p> FixedVolumeConductionLink<'l, 'p> {
    /// Bind the existing ordered interface without replacing any timeline state.
    pub fn try_aligned(
        link: &ThermalLink,
        laws: [&'l FixedVolumeThermalLaw<'p>; 2],
        distances_m: [f64; 2],
        contact_resistance_m2_k_w: f64,
    ) -> Result<Self, DepositionError> {
        if link.a == link.b
            || !positive(distances_m[0])
            || !positive(distances_m[1])
            || !nonnegative(link.area_m2)
            || !nonnegative(contact_resistance_m2_k_w)
        {
            return Err(DepositionError::InvalidInput(
                "invalid sourced conduction interface geometry",
            ));
        }
        Ok(Self {
            id: link.id,
            cells: [link.a, link.b],
            laws,
            distances_m,
            area_m2: link.area_m2,
            contact_resistance_m2_k_w,
        })
    }

    #[must_use]
    pub const fn id(&self) -> u64 {
        self.id
    }

    /// Evaluate only the matching original laws, cell bases and interface geometry.
    /// Return this value in the existing `StepCoefficients::conductances_w_k`
    /// during an implicit endpoint solve, or explicitly refresh an existing link.
    /// Failed evaluation publishes no state and never freezes a missing source.
    pub fn conductance_w_k(
        &self,
        state: &DepositionSimulation<'_>,
    ) -> Result<f64, DepositionError> {
        let link = state
            .links()
            .get(&self.id)
            .ok_or(DepositionError::InvalidInput(
                "sourced conduction link is not active",
            ))?;
        if [link.a, link.b] != self.cells || link.area_m2 != self.area_m2 {
            return Err(DepositionError::InvalidInput(
                "sourced conduction interface geometry changed",
            ));
        }
        let mut conductivity = [0.0; 2];
        for (i, value) in conductivity.iter_mut().enumerate() {
            let cell = state
                .cells()
                .get(&self.cells[i])
                .ok_or(DepositionError::InvalidInput(
                    "sourced conduction cell is not active",
                ))?;
            if !state.uses_material_law(cell.material_id, self.laws[i]) {
                return Err(DepositionError::InvalidInput(
                    "sourced conduction requires the original caloric law",
                ));
            }
            self.laws[i].admit_cell(cell)?;
            let temperature = state.temperature_k(cell.id)?;
            *value = self.laws[i]
                .path
                .at_temperature(temperature)
                .map_err(|_| DepositionError::MaterialDomain)?
                .conductivity_w_m_k();
        }
        series_contact_conductance_w_k(
            self.area_m2,
            self.distances_m,
            conductivity,
            self.contact_resistance_m2_k_w,
        )
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn invalid(message: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(message)
}
