//! Bounded equilibrium enthalpy/porosity closure for a stationary drag frame.
//!
//! The permeability scale and admitted fraction are explicit model data, not
//! universal consequences of melting. This single-velocity mixture closure
//! supplies neither solid elasticity nor separate constituent velocities.
//! The caller must justify the chart and closure for one fixed pressure and
//! composition context. Raw model data acquire no measured source authority.

use fs_blake3::{ContentHash, DomainHasher};

use super::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use crate::MaterialError;

/// Phase volumes per unit mixture mass, converted with constituent densities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseVolumeFractions {
    pub solid: f64,
    pub liquid: f64,
    pub bulk_density_kg_m3: f64,
}

impl PhaseVolumeFractions {
    /// No equality of mass and volume fractions is presumed. The supplied
    /// phase densities must also reproduce the original chart's bulk density.
    pub fn try_new(
        state: EquilibriumPhaseState,
        solid_density_kg_m3: f64,
        liquid_density_kg_m3: f64,
    ) -> Result<Self, MaterialError> {
        if ![solid_density_kg_m3, liquid_density_kg_m3]
            .into_iter()
            .all(|rho| rho.is_finite() && rho > 0.0)
        {
            return Err(invalid("constituent densities must be positive and finite"));
        }
        let solid = state.solid_mass_fraction() / solid_density_kg_m3;
        let liquid = state.liquid_mass_fraction() / liquid_density_kg_m3;
        let volume = solid + liquid;
        let rho = if solid_density_kg_m3.to_bits() == liquid_density_kg_m3.to_bits() {
            solid_density_kg_m3
        } else {
            volume.recip()
        };
        if !volume.is_finite() || volume <= 0.0 || !rho.is_finite() || rho <= 0.0 {
            return Err(invalid("mixture specific volume is not representable"));
        }
        if (rho - state.bulk_density_kg_m3()).abs()
            > 64.0 * f64::EPSILON * rho.max(state.bulk_density_kg_m3())
        {
            return Err(invalid(
                "constituent volumes disagree with the original enthalpy chart density",
            ));
        }
        let liquid = if solid_density_kg_m3.to_bits() == liquid_density_kg_m3.to_bits() {
            state.liquid_mass_fraction()
        } else if state.liquid_mass_fraction() == 0.0 {
            0.0
        } else if state.liquid_mass_fraction() == 1.0 {
            1.0
        } else {
            liquid / volume
        };
        if !liquid.is_finite()
            || !(0.0..=1.0).contains(&liquid)
            || (state.liquid_mass_fraction() > 0.0 && liquid == 0.0)
            || (state.solid_mass_fraction() > 0.0 && liquid == 1.0)
        {
            return Err(invalid("mixture volume fraction is not representable"));
        }
        Ok(Self {
            solid: 1.0 - liquid,
            liquid,
            bulk_density_kg_m3: rho,
        })
    }
}

/// Explicit Carman–Kozeny closure `K^-1=(1-phi)^2/(K0*phi^3)`.
///
/// `K0` is a positive model length squared, not a fitted universal constant.
/// No artificial epsilon makes fully solid material weakly mobile. The
/// original constructor admits only its positive declared liquid-fraction
/// range. The explicit immobile-solid constructor adds an exact stationary
/// solid branch and positive fractions with representable finite mobility. Fully liquid
/// material has zero extra drag. The original liquid viscosity later gives
/// `drag=mu/K` in N s/m4. This stationary, homogeneous, single-velocity
/// enthalpy/porosity approximation does not provide a moving elastic skeleton,
/// two-velocity phase transport or a stress-free resolidification reference.
#[derive(Clone, Debug)]
pub struct MushyPhaseLaw<'a> {
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    solid_density_kg_m3: f64,
    liquid_density_kg_m3: f64,
    permeability_scale_m2: f64,
    minimum_liquid_volume_fraction: f64,
    identity: ContentHash,
}

/// Same accepted phase state with its declared volume and mobility response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MushyPhaseState {
    pub phase: EquilibriumPhaseState,
    pub volumes: PhaseVolumeFractions,
    /// Inverse permeability [1/m2], zero at the fully liquid endpoint.
    pub inverse_permeability_per_m2: f64,
    pub law_identity: ContentHash,
}

/// Exact stationary solid constraint or a finite mobile phase response.
///
/// Immobility is a separate constitutive branch, never an infinite viscosity,
/// clipped fraction or large finite drag. A momentum owner must enforce its
/// zero velocity while keeping the original thermal mass and source state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PhaseMobility {
    ImmobileSolid {
        phase: EquilibriumPhaseState,
        volumes: PhaseVolumeFractions,
        law_identity: ContentHash,
    },
    Mobile(MushyPhaseState),
}
impl PhaseMobility {
    #[must_use]
    pub const fn phase(&self) -> EquilibriumPhaseState {
        match self {
            Self::ImmobileSolid { phase, .. } => *phase,
            Self::Mobile(state) => state.phase,
        }
    }
    #[must_use]
    pub const fn volumes(&self) -> PhaseVolumeFractions {
        match self {
            Self::ImmobileSolid { volumes, .. } => *volumes,
            Self::Mobile(state) => state.volumes,
        }
    }
    #[must_use]
    pub const fn law_identity(&self) -> ContentHash {
        match self {
            Self::ImmobileSolid { law_identity, .. } => *law_identity,
            Self::Mobile(state) => state.law_identity,
        }
    }
    #[must_use]
    pub const fn mobile(&self) -> Option<&MushyPhaseState> {
        match self {
            Self::ImmobileSolid { .. } => None,
            Self::Mobile(state) => Some(state),
        }
    }
}

impl<'a> MushyPhaseLaw<'a> {
    /// Bind the exact chart and constant phase densities. Checking every knot
    /// suffices for its piecewise-linear specific-volume/fraction convention.
    /// The entire curve may include solid states; this mobility model admits
    /// only states meeting its separately declared liquid-fraction minimum.
    pub fn try_new(
        curve: &'a EquilibriumEnthalpyPhaseCurve,
        solid_density_kg_m3: f64,
        liquid_density_kg_m3: f64,
        permeability_scale_m2: f64,
        minimum_liquid_volume_fraction: f64,
    ) -> Result<Self, MaterialError> {
        if !permeability_scale_m2.is_finite()
            || permeability_scale_m2 <= 0.0
            || !minimum_liquid_volume_fraction.is_finite()
            || minimum_liquid_volume_fraction <= 0.0
            || minimum_liquid_volume_fraction >= 1.0
        {
            return Err(invalid(
                "mushy permeability scale must be positive and minimum liquid volume fraction in (0,1)",
            ));
        }
        Self::build(
            curve,
            solid_density_kg_m3,
            liquid_density_kg_m3,
            permeability_scale_m2,
            minimum_liquid_volume_fraction,
        )
    }

    /// Admit the exact fully solid endpoint as a stationary constraint and all
    /// strictly positive liquid fractions whose Carman–Kozeny resistance is
    /// representable. There is no positive cutoff, epsilon, or fraction clipping.
    ///
    /// This supplies no elastic skeleton, residual stress or shrinkage law.
    /// The consuming momentum owner must implement the exact immobile branch
    /// and preserve thermal mass/enthalpy when entering or leaving it.
    pub fn try_new_with_immobile_solid(
        curve: &'a EquilibriumEnthalpyPhaseCurve,
        solid_density_kg_m3: f64,
        liquid_density_kg_m3: f64,
        permeability_scale_m2: f64,
    ) -> Result<Self, MaterialError> {
        if !permeability_scale_m2.is_finite() || permeability_scale_m2 <= 0.0 {
            return Err(invalid(
                "mushy permeability scale must be positive and finite",
            ));
        }
        Self::build(
            curve,
            solid_density_kg_m3,
            liquid_density_kg_m3,
            permeability_scale_m2,
            0.0,
        )
    }

    fn build(
        curve: &'a EquilibriumEnthalpyPhaseCurve,
        solid_density_kg_m3: f64,
        liquid_density_kg_m3: f64,
        permeability_scale_m2: f64,
        minimum_liquid_volume_fraction: f64,
    ) -> Result<Self, MaterialError> {
        for knot in curve.knots() {
            let state = curve
                .state_at_specific_enthalpy(knot.specific_enthalpy_j_kg)
                .map_err(|error| MaterialError::State {
                    what: error.to_string(),
                })?;
            PhaseVolumeFractions::try_new(state, solid_density_kg_m3, liquid_density_kg_m3)?;
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.mushy-carman-kozeny.v1");
        hash.update(curve.identity().as_bytes());
        for value in [
            solid_density_kg_m3,
            liquid_density_kg_m3,
            permeability_scale_m2,
            minimum_liquid_volume_fraction,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        let law = Self {
            curve,
            solid_density_kg_m3,
            liquid_density_kg_m3,
            permeability_scale_m2,
            minimum_liquid_volume_fraction,
            identity: hash.finalize(),
        };
        if minimum_liquid_volume_fraction > 0.0 {
            law.inverse_permeability(minimum_liquid_volume_fraction)?;
        }
        Ok(law)
    }

    pub const fn curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.curve
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Whether exact zero liquid fraction has an explicit stationary branch.
    #[must_use]
    pub const fn admits_immobile_solid(&self) -> bool {
        self.minimum_liquid_volume_fraction == 0.0
    }

    /// Only this exact constant-density subset can enter a fixed-volume
    /// incompressible transport owner without a phase-volume source.
    pub fn fixed_density_kg_m3(&self) -> Option<f64> {
        let rho = self.solid_density_kg_m3;
        (rho.to_bits() == self.liquid_density_kg_m3.to_bits()
            && self
                .curve
                .knots()
                .iter()
                .all(|k| k.bulk_density_kg_m3.to_bits() == rho.to_bits()))
        .then_some(rho)
    }

    /// Resolve from h, never temperature alone, so latent-plateau melting
    /// changes mobility without inventing a heat capacity or adding latent heat.
    pub fn at_specific_enthalpy(&self, h: f64) -> Result<MushyPhaseState, MaterialError> {
        match self.classify_at_specific_enthalpy(h)? {
            PhaseMobility::Mobile(state) => Ok(state),
            PhaseMobility::ImmobileSolid { .. } => Err(invalid(
                "exactly solid material has no finite Darcy response; use the immobile phase branch",
            )),
        }
    }

    /// Resolve topology before any liquid-only property query. The original
    /// positive-cutoff constructor retains its admitted fraction domain; only
    /// the explicit immobile-solid constructor admits the exact solid branch.
    /// Finite positive fractions never become solid through a numerical cutoff.
    pub fn classify_at_specific_enthalpy(&self, h: f64) -> Result<PhaseMobility, MaterialError> {
        let phase =
            self.curve
                .state_at_specific_enthalpy(h)
                .map_err(|error| MaterialError::State {
                    what: error.to_string(),
                })?;
        let volumes = PhaseVolumeFractions::try_new(
            phase,
            self.solid_density_kg_m3,
            self.liquid_density_kg_m3,
        )?;
        if volumes.liquid == 0.0 && self.minimum_liquid_volume_fraction == 0.0 {
            return Ok(PhaseMobility::ImmobileSolid {
                phase,
                volumes,
                law_identity: self.identity,
            });
        }
        if volumes.liquid < self.minimum_liquid_volume_fraction {
            return Err(invalid(
                "mushy state is below its admitted liquid volume fraction; solid locking needs its own owner",
            ));
        }
        Ok(PhaseMobility::Mobile(MushyPhaseState {
            phase,
            volumes,
            inverse_permeability_per_m2: self.inverse_permeability(volumes.liquid)?,
            law_identity: self.identity,
        }))
    }

    fn inverse_permeability(&self, phi: f64) -> Result<f64, MaterialError> {
        if phi == 1.0 {
            return Ok(0.0);
        }
        let ratio = (1.0 - phi) / phi;
        let direct = (ratio * ratio / phi) / self.permeability_scale_m2;
        let inverse = if direct.is_finite() && direct > 0.0 {
            direct
        } else {
            // Preserve the established arithmetic in ordinary cases. Scaling
            // first avoids intermediate overflow when a large K0 compensates
            // a small positive fraction. An unrepresentable final resistance
            // still refuses; no fraction floor or drag cap is introduced.
            let scaled = (1.0 - phi) / (fs_math::det::sqrt(self.permeability_scale_m2) * phi);
            (scaled * scaled) / phi
        };
        if !inverse.is_finite() || inverse <= 0.0 {
            return Err(invalid("mushy inverse permeability is not representable"));
        }
        Ok(inverse)
    }
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::State { what: what.into() }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::phase::EnthalpyPhaseKnot;

    fn curve(rhos: [f64; 2]) -> EquilibriumEnthalpyPhaseCurve {
        EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([8; 32]),
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: 500.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: rhos[0],
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 100.0,
                    temperature_k: 500.0,
                    liquid_mass_fraction: 1.0,
                    bulk_density_kg_m3: rhos[1],
                },
            ],
        )
        .unwrap()
    }

    #[test]
    fn unequal_density_converts_mass_to_volume_and_retains_latent_state() {
        let curve = curve([2.0, 1.0]);
        let law = MushyPhaseLaw::try_new(&curve, 2.0, 1.0, 2.0, 0.01).unwrap();
        let state = law.at_specific_enthalpy(50.0).unwrap();
        assert_eq!(state.phase.liquid_mass_fraction(), 0.5);
        assert!((state.volumes.liquid - 2.0 / 3.0).abs() < 1e-15);
        assert!((state.volumes.bulk_density_kg_m3 - 4.0 / 3.0).abs() < 1e-15);
        assert!((state.inverse_permeability_per_m2 - 3.0 / 16.0).abs() < 1e-15);
        assert_eq!(state.phase.temperature_k(), 500.0);
        assert_eq!(state.phase.specific_enthalpy_j_kg(), 50.0);
        assert_eq!(law.fixed_density_kg_m3(), None);
        assert!(law.at_specific_enthalpy(0.0).is_err());
        let liquid = law.at_specific_enthalpy(100.0).unwrap();
        assert_eq!(liquid.volumes.liquid, 1.0);
        assert_eq!(liquid.inverse_permeability_per_m2, 0.0);
        let solid =
            PhaseVolumeFractions::try_new(curve.state_at_specific_enthalpy(0.0).unwrap(), 2.0, 1.0)
                .unwrap();
        assert_eq!(solid.solid, 1.0);
    }

    #[test]
    fn density_and_mobility_support_are_explicit_and_no_solid_leak_is_invented() {
        let curve = curve([2.0, 2.0]);
        assert!(MushyPhaseLaw::try_new(&curve, 2.0, 1.0, 1.0, 0.1).is_err());
        assert!(MushyPhaseLaw::try_new(&curve, 2.0, 2.0, 0.0, 0.1).is_err());
        let law = MushyPhaseLaw::try_new(&curve, 2.0, 2.0, 1.0, 0.1).unwrap();
        assert_eq!(law.fixed_density_kg_m3(), Some(2.0));
        assert!(law.at_specific_enthalpy(9.0).is_err());
        assert!(law.at_specific_enthalpy(101.0).is_err());
        assert!(law.at_specific_enthalpy(f64::NAN).is_err());
        assert_eq!(
            law.at_specific_enthalpy(50.0)
                .unwrap()
                .inverse_permeability_per_m2,
            2.0
        );
        assert!(
            law.at_specific_enthalpy(75.0)
                .unwrap()
                .inverse_permeability_per_m2
                < 2.0
        );
    }
    #[test]
    fn explicit_solid_branch_has_no_positive_cutoff_or_finite_drag_surrogate() {
        let curve = curve([2.0, 2.0]);
        let strict = MushyPhaseLaw::try_new(&curve, 2.0, 2.0, 1.0, 0.1).unwrap();
        let full = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 1.0).unwrap();
        assert!(!strict.admits_immobile_solid());
        assert!(full.admits_immobile_solid());
        assert_ne!(strict.identity(), full.identity());
        assert!(strict.classify_at_specific_enthalpy(0.0).is_err());
        assert!(strict.classify_at_specific_enthalpy(9.0).is_err());
        assert!(MushyPhaseLaw::try_new(&curve, 2.0, 2.0, 1.0, 0.0).is_err());
        let solid = full.classify_at_specific_enthalpy(0.0).unwrap();
        assert!(matches!(solid, PhaseMobility::ImmobileSolid { .. }));
        assert_eq!(
            solid.phase(),
            curve.state_at_specific_enthalpy(0.0).unwrap()
        );
        assert_eq!(solid.volumes().liquid, 0.0);
        assert_eq!(solid.mobile(), None);
        // The legacy API must not smuggle solid topology into a Darcy tuple.
        assert!(full.at_specific_enthalpy(0.0).is_err());
        for h in [1e-10, 9.0, 50.0, 100.0] {
            let state = full.classify_at_specific_enthalpy(h).unwrap();
            assert!(state.mobile().is_some());
            assert!(state.phase().liquid_mass_fraction() > 0.0);
        }
        assert_eq!(
            full.at_specific_enthalpy(100.0)
                .unwrap()
                .inverse_permeability_per_m2,
            0.0
        );
    }

    #[test]
    fn compensating_permeability_scaling_is_finite_but_overflow_never_becomes_solid() {
        let curve = curve([2.0, 2.0]);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 1e300).unwrap();
        // The direct phi^-3 intermediate overflows, although the final value
        // after the explicit model K0 is finite.
        let state = law.at_specific_enthalpy(1e-108).unwrap();
        assert!((state.inverse_permeability_per_m2 / 1e30 - 1.0).abs() < 1e-14);
        let ordinary = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 2.0, 2.0, 1.0).unwrap();
        assert!(ordinary.classify_at_specific_enthalpy(1e-108).is_err());
        assert!(
            ordinary
                .classify_at_specific_enthalpy(0.0)
                .unwrap()
                .mobile()
                .is_none()
        );
    }
}
