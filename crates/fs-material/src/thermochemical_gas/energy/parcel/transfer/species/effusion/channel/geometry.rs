//! Geometry-only transmission for a straight circular molecular channel.
//!
//! Santeler's approximation is an **Estimated** Clausing transmission, not an
//! exact transport solution or a certified error enclosure. With `h = L/D`,
//! it is `(14 + 4h)/(14 + 18h + 3h²)`. See Chiggiato, *Vacuum Technology for
//! Ion Sources*, Eq. (21), <https://arxiv.org/abs/1404.0960>, and Santeler,
//! <https://doi.org/10.1116/1.573923>. The direct probability is the exact
//! coaxial equal-disc view factor `(sqrt(1+h²) - h)²` for cosine incidence.
//!
//! The model assumes quiescent Maxwellian reservoirs, nonadsorbing diffuse
//! walls, and collisionless trajectories. The thermal consumer additionally
//! owns a stationary, isothermal wall and qualifies its energy accommodation.
//! Neither arbitrary incidence/drift nor specular reflection can reuse these
//! probabilities.

mod clausing;
pub use clausing::accommodation::{AccommodationResolution, MolecularAccommodation};
pub use clausing::momentum::{
    AxialAccommodationResolution, AxialMomentumResolution, MolecularAxialAccommodation,
    MolecularAxialMoments,
};
pub use clausing::{ClausingControls, ClausingResolution};

use super::super::{GasEnergyError, invalid, positive};
use fs_math::det;

/// A qualified circular bore, including a stable partition of molecular paths.
///
/// Diameter and length must each be at most one tenth of the authored minimum
/// mean free path. The bound must cover all participating gas mixtures and
/// temperatures/pressures, including the wall temperature. These are admission
/// gates, not a measured accuracy guarantee. Total exposed area is supplied
/// separately by the consumer, permitting multiple identical bores.
#[derive(Clone, Debug, PartialEq)]
pub struct CylindricalMolecularGeometry {
    diameter_m: f64,
    length_m: f64,
    minimum_mean_free_path_m: f64,
    provenance: String,
    transmission: f64,
    direct: f64,
    wall_hit: f64,
    reflected: f64,
    wall_transmitted: f64,
    resolution: Option<ClausingResolution>,
}

impl CylindricalMolecularGeometry {
    pub fn try_new(
        diameter_m: f64,
        length_m: f64,
        minimum_mean_free_path_m: f64,
        provenance: impl Into<String>,
    ) -> Result<Self, GasEnergyError> {
        let provenance = provenance.into();
        if !positive(diameter_m)
            || !positive(minimum_mean_free_path_m)
            || !length_m.is_finite()
            || length_m < 0.0
            || diameter_m / minimum_mean_free_path_m > 0.1
            || length_m / minimum_mean_free_path_m > 0.1
            || provenance.trim().is_empty()
            || provenance.len() > 4096
        {
            return Err(invalid(
                "molecular channel needs finite dimensions, collisionless diameter and length, and bounded qualification provenance",
            ));
        }
        let (transmission, direct, wall_hit, reflected, wall_transmitted) = if length_m == 0.0 {
            (1.0, 1.0, 0.0, 0.0, 0.0)
        } else {
            let h = length_m / diameter_m;
            if !positive(h) {
                return Err(invalid(
                    "molecular channel aspect ratio is not representable",
                ));
            }
            let (transmission, direct, wall_hit, reflected, wall_transmitted) = if h <= 1.0 {
                // Retain the small path fractions directly: subtracting two
                // probabilities rounded to one would erase finite wall heat.
                let wall_hit = 2.0 * h / (det::sqrt(1.0 + h * h) + h);
                let reflected = h * ((14.0 + 3.0 * h) / (14.0 + h * (18.0 + 3.0 * h)));
                (
                    1.0 - reflected,
                    1.0 - wall_hit,
                    wall_hit,
                    reflected,
                    wall_hit - reflected,
                )
            } else {
                // Reciprocal aspect ratio avoids h² overflow and preserves
                // long-channel direct and transmitted paths near zero.
                let inverse_h = 1.0 / h;
                let transmission = inverse_h * (4.0 + 14.0 * inverse_h)
                    / (3.0 + inverse_h * (18.0 + 14.0 * inverse_h));
                let direct_root = inverse_h / (det::sqrt(1.0 + inverse_h * inverse_h) + 1.0);
                let direct = direct_root * direct_root;
                (
                    transmission,
                    direct,
                    1.0 - direct,
                    1.0 - transmission,
                    transmission - direct,
                )
            };
            if ![transmission, direct, wall_hit, reflected, wall_transmitted]
                .iter()
                .all(|&probability| positive(probability) && probability <= 1.0)
                || direct > transmission
            {
                return Err(invalid(
                    "molecular channel path probabilities are not representable",
                ));
            }
            (transmission, direct, wall_hit, reflected, wall_transmitted)
        };
        Ok(Self {
            diameter_m,
            length_m,
            minimum_mean_free_path_m,
            provenance,
            transmission,
            direct,
            wall_hit,
            reflected,
            wall_transmitted,
            resolution: None,
        })
    }

    /// Resolve the diffuse-wall Clausing equation on successively refined axial
    /// rings, retaining the original geometry and collision qualification.
    ///
    /// The returned probabilities feed every existing channel consumer. The
    /// refinement estimate is numerical evidence, not a certified enclosure or
    /// an extension of the physical diffuse/collisionless applicability domain.
    pub fn resolve_transmission(
        self,
        controls: ClausingControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasEnergyError> {
        clausing::resolve(self, controls, cancelled)
    }

    /// Compose constant per-encounter energy accommodation with all repeated
    /// diffuse wall collisions. The coefficient does not change scattering
    /// angles or mass transport. Interior coefficients refine the geometry and
    /// both energy weights together; exact endpoint laws retain this geometry.
    ///
    /// A material consumer owns qualification of this coefficient over its
    /// original gas catalog and caloric support and binds the returned pair.
    /// The wall is stationary and isothermal, with no adsorption or reactions.
    pub fn resolve_thermal_accommodation(
        self,
        coefficient: f64,
        controls: ClausingControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, MolecularAccommodation), GasEnergyError> {
        clausing::accommodation::resolve(self, coefficient, controls, cancelled)
    }

    /// Resolve one to `MAX_EFFUSIVE_SPECIES` energy-accommodation coefficients
    /// on a shared Clausing grid and mass solution. Results preserve input order;
    /// identical coefficients share thermal solves. Refinement must satisfy both
    /// heat weights for every coefficient and the common mass-transmission law.
    /// Pure endpoint batches preserve the supplied geometry. In a mixed batch,
    /// coefficient one uses the final refined geometry's wall-path probabilities.
    pub fn resolve_thermal_accommodation_batch(
        self,
        coefficients: &[f64],
        controls: ClausingControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, Vec<MolecularAccommodation>), GasEnergyError> {
        clausing::accommodation::resolve_batch(self, coefficients, controls, cancelled)
    }

    /// Resolve the direct, wall-transmitted and returned axial angular moments
    /// together with mass transmission on one Clausing grid. Each moment is the
    /// expectation of the nonnegative outgoing axial direction cosine times
    /// its path indicator, normalized by the ORIGINAL incident molecule count.
    ///
    /// These are geometric weights, not speeds or forces. A material consumer
    /// must separately qualify the incident and reflected speed distributions;
    /// an energy-accommodation coefficient alone cannot supply that law.
    pub fn resolve_axial_momentum(
        self,
        controls: ClausingControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, MolecularAxialMoments), GasEnergyError> {
        clausing::momentum::resolve(self, controls, cancelled)
    }

    /// Resolve mass, heat and angular/collision-history moments on the same
    /// final grid for one to `MAX_EFFUSIVE_SPECIES` coefficients. Each wall
    /// encounter independently resets the molecular speed and internal state
    /// to the isothermal wall distribution with probability alpha; otherwise
    /// it retains them. Both outcomes scatter diffusely and independently of
    /// speed. This explicit Bernoulli law is stronger than mean energy
    /// accommodation alone and must be qualified by the material consumer.
    ///
    /// Results retain input order and duplicate coefficients reuse their solve.
    /// Retained and reset angular weights partition the total outgoing moments;
    /// each response also owns its heat probabilities on these exact rings.
    /// One mass solve and one shared angular-escape pair per grid serve every
    /// distinct interior coefficient, with O(rings + species) scratch storage.
    pub fn resolve_axial_accommodation_batch(
        self,
        coefficients: &[f64],
        controls: ClausingControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        (
            Self,
            MolecularAxialMoments,
            Vec<MolecularAxialAccommodation>,
        ),
        GasEnergyError,
    > {
        clausing::momentum::resolve_accommodation_batch(self, coefficients, controls, cancelled)
    }

    /// `None` identifies the original Santeler approximation. A retained report
    /// identifies the actual ring refinement used for these probabilities.
    pub const fn resolution(&self) -> Option<&ClausingResolution> {
        self.resolution.as_ref()
    }

    pub const fn diameter_m(&self) -> f64 {
        self.diameter_m
    }
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    pub const fn minimum_mean_free_path_m(&self) -> f64 {
        self.minimum_mean_free_path_m
    }
    pub fn provenance(&self) -> &str {
        &self.provenance
    }
    /// Reciprocal transmission: Santeler by default, or the explicitly resolved
    /// Clausing value when `resolution()` is present.
    pub const fn transmission_probability(&self) -> f64 {
        self.transmission
    }
    /// Exact geometric view factor, evaluated to floating-point precision.
    pub const fn direct_probability(&self) -> f64 {
        self.direct
    }
    /// Probability of at least one wall collision, retaining tiny bore lengths.
    pub const fn wall_hit_probability(&self) -> f64 {
        self.wall_hit
    }
    /// Probability of returning to the entrance under the selected geometry law.
    pub const fn reflected_probability(&self) -> f64 {
        self.reflected
    }
    /// Probability of transmission after wall impact under the selected law.
    pub const fn wall_transmitted_probability(&self) -> f64 {
        self.wall_transmitted
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn geometry(h: f64) -> CylindricalMolecularGeometry {
        CylindricalMolecularGeometry::try_new(1.0, h, 100.0 * h.max(1.0), "test qualification")
            .unwrap()
    }

    #[test]
    fn zero_length_is_the_thin_aperture_limit() {
        let g = geometry(0.0);
        assert_eq!(g.transmission_probability(), 1.0);
        assert_eq!(g.direct_probability(), 1.0);
        assert_eq!(g.wall_hit_probability(), 0.0);
        assert_eq!(g.reflected_probability(), 0.0);
        assert_eq!(g.wall_transmitted_probability(), 0.0);
    }

    #[test]
    fn santeler_agrees_with_independent_cole_values_within_its_estimate() {
        // Cole's L/R table, reproduced in COMSOL's vacuum_capillary benchmark:
        // https://doc.comsol.com/6.4/doc/com.comsol.help.models.molec.vacuum_capillary/vacuum_capillary.html
        for (h, cole) in [
            (1.0, 0.5142),
            (2.0, 0.3566),
            (3.0, 0.2755),
            (4.0, 0.2253),
            (5.0, 0.1910),
        ] {
            let g = geometry(h);
            assert!((g.transmission_probability() / cole - 1.0).abs() < 0.007);
            assert!(
                (g.direct_probability()
                    + g.wall_transmitted_probability()
                    + g.reflected_probability()
                    - 1.0)
                    .abs()
                    < 4.0 * f64::EPSILON
            );
        }
        // The equal-diameter/length view factor is 3 - 2 sqrt(2).
        assert!((geometry(1.0).direct_probability() - 0.171_572_875_253_809_9).abs() < 1e-15);
    }

    #[test]
    fn short_channels_retain_wall_heat_when_direct_probability_rounds_to_one() {
        let h = 1e-18;
        let g = geometry(h);
        assert_eq!(g.direct_probability(), 1.0);
        assert_eq!(g.transmission_probability(), 1.0);
        assert!((g.wall_hit_probability() / h - 2.0).abs() < 1e-15);
        assert!((g.reflected_probability() / h - 1.0).abs() < 1e-15);
        assert!((g.wall_transmitted_probability() / h - 1.0).abs() < 1e-15);
    }

    #[test]
    fn long_channels_preserve_knudsen_and_direct_view_factor_asymptotes() {
        for h in [1e8, 1e150] {
            let g = geometry(h);
            assert!((g.transmission_probability() * h - 4.0 / 3.0).abs() < 1e-7);
            assert!((g.direct_probability() * h * h - 0.25).abs() < 1e-15);
            assert!(g.wall_transmitted_probability() > 0.0);
        }
    }

    #[test]
    fn rejects_invalid_qualification_and_unrepresentable_paths() {
        for (d, l, lambda) in [
            (0.0, 1.0, 100.0),
            (1.0, -1.0, 100.0),
            (1.0, f64::NAN, 100.0),
            (f64::INFINITY, 1.0, 100.0),
            (1.0, 1.0, f64::INFINITY),
            (1.0, 0.0, 1.0),
            (0.001, 1.0, 1.0),
            (1e200, 1e-200, 1e202),
            (1e-200, 1.0, 100.0),
        ] {
            assert!(CylindricalMolecularGeometry::try_new(d, l, lambda, "qualification").is_err());
        }
        assert!(CylindricalMolecularGeometry::try_new(1.0, 1.0, 100.0, " ").is_err());
        assert!(CylindricalMolecularGeometry::try_new(1.0, 1.0, 100.0, "x".repeat(4097)).is_err());
    }
}
