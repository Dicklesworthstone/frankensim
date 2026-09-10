//! Equilibrium solid-liquid state resolved on a specific-enthalpy coordinate.
//!
//! Total enthalpy is the primary coordinate because it remains single valued
//! through an isothermal phase-change plateau. A caller supplies a monotone,
//! evidence-derived curve; this module neither recognizes material names nor
//! invents melting points, latent heats, densities, or phase fractions.
//!
//! This is a constitutive state primitive, not a heat-transfer or free-surface
//! solver. Thermal transport owns changes in specific enthalpy. Solid,
//! finite-strain, fluid, remeshing, acoustic, and optical consumers own their
//! respective responses to the returned phase state.

use core::fmt;

use fs_blake3::{ContentHash, DomainHasher};

const EQUILIBRIUM_ENTHALPY_PHASE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.equilibrium-enthalpy-phase.v1";
const EQUILIBRIUM_PHASE_STATE_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-material.equilibrium-phase-state.v1";

/// Proposed uniform-body thermal step shared by transport and coupled owners.
/// Geometry and mass belong to the specimen; this carrier does not justify
/// lumping, solve transport, or mutate the accepted state.
#[derive(Clone, Copy, Debug)]
pub struct UniformEnthalpyStepInput<'a> {
    /// Immutable equilibrium chart used by the accepted body.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    /// Accepted initial state on that chart.
    pub initial: EquilibriumPhaseState,
    /// Invariant specimen mass [kg].
    pub mass_kg: f64,
    /// Specimen volume [m3].
    pub volume_m3: f64,
    /// Whole exposed boundary area [m2].
    pub surface_area_m2: f64,
    /// Internally deposited heat over this step [J], counted exactly once.
    pub internal_heat_j: f64,
    /// Physical duration of this step [s].
    pub duration_s: f64,
}

/// Proposed thermal state and boundary transfer, with a solver-owned report.
/// The coupled owner must check the chart, energy balance and its own regime
/// before publishing this proposal. A transport callback must not debit a
/// mutable external reservoir before the coupled transaction accepts it.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformEnthalpyStep<R> {
    /// Proposed state on the input equilibrium chart.
    pub state: EquilibriumPhaseState,
    /// Signed integrated boundary heat into the body [J].
    pub external_heat_j: f64,
    /// Declared transport solve residual allowance [J], separate from roundoff
    /// and from temporal or constitutive-model error.
    pub energy_residual_tolerance_j: f64,
    /// Actual transport result retained without erasing its concrete type.
    pub report: R,
}

/// One source-provided point on a solid-liquid equilibrium enthalpy curve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnthalpyPhaseKnot {
    /// Specific enthalpy [J/kg]. Knots must be strictly increasing in this
    /// coordinate.
    pub specific_enthalpy_j_kg: f64,
    /// Absolute temperature [K]. It must be positive and nondecreasing; equal
    /// temperatures are deliberately admitted for an isothermal latent-heat
    /// plateau.
    pub temperature_k: f64,
    /// Equilibrium liquid mass fraction in `[0, 1]`, nondecreasing with
    /// enthalpy. The solid mass fraction is exactly its complement.
    pub liquid_mass_fraction: f64,
    /// Equilibrium bulk density [kg/m3] at this state.
    pub bulk_density_kg_m3: f64,
}

/// Coarse phase topology selected from an equilibrium mass fraction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolidLiquidPhase {
    /// Exactly zero liquid mass fraction.
    Solid,
    /// Coexisting solid and liquid mass fractions.
    SolidLiquid,
    /// Exactly unit liquid mass fraction.
    Liquid,
}

/// State resolved from one admitted specific enthalpy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumPhaseState {
    /// Specific enthalpy [J/kg].
    specific_enthalpy_j_kg: f64,
    /// Absolute temperature [K].
    temperature_k: f64,
    /// Solid mass fraction in `[0, 1]`.
    solid_mass_fraction: f64,
    /// Liquid mass fraction in `[0, 1]`.
    liquid_mass_fraction: f64,
    /// Equilibrium bulk density [kg/m3].
    bulk_density_kg_m3: f64,
    /// Exact coarse phase topology.
    phase: SolidLiquidPhase,
    /// Material card that owns the equilibrium curve.
    material_card_identity: ContentHash,
    /// Identity of the complete admitted equilibrium curve.
    phase_curve_identity: ContentHash,
    /// Identity binding the curve and this resolved state.
    identity: ContentHash,
}

impl EquilibriumPhaseState {
    /// Specific enthalpy [J/kg].
    #[must_use]
    pub const fn specific_enthalpy_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg
    }

    /// Absolute temperature [K].
    #[must_use]
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }

    /// Solid mass fraction in `[0, 1]`.
    #[must_use]
    pub const fn solid_mass_fraction(self) -> f64 {
        self.solid_mass_fraction
    }

    /// Liquid mass fraction in `[0, 1]`.
    #[must_use]
    pub const fn liquid_mass_fraction(self) -> f64 {
        self.liquid_mass_fraction
    }

    /// Equilibrium bulk density [kg/m3].
    #[must_use]
    pub const fn bulk_density_kg_m3(self) -> f64 {
        self.bulk_density_kg_m3
    }

    /// Exact coarse phase topology.
    #[must_use]
    pub const fn phase(self) -> SolidLiquidPhase {
        self.phase
    }

    /// Material card that owns the equilibrium curve.
    #[must_use]
    pub const fn material_card_identity(self) -> ContentHash {
        self.material_card_identity
    }

    /// Identity of the complete admitted equilibrium curve.
    #[must_use]
    pub const fn phase_curve_identity(self) -> ContentHash {
        self.phase_curve_identity
    }

    /// Identity binding the curve and this resolved state.
    #[must_use]
    pub const fn identity(self) -> ContentHash {
        self.identity
    }

    /// Whether a solid-only constitutive/dynamics rung may consider this state.
    ///
    /// This is necessary but not sufficient admission: a solid consumer must
    /// still resolve its own temperature-, pressure-, rate-, and history-
    /// dependent properties at the same state.
    #[must_use]
    pub const fn is_fully_solid(self) -> bool {
        matches!(self.phase, SolidLiquidPhase::Solid)
    }
}

/// Immutable equilibrium phase curve bound to one material-card identity.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumEnthalpyPhaseCurve {
    material_card_identity: ContentHash,
    knots: Vec<EnthalpyPhaseKnot>,
    identity: ContentHash,
}

impl EquilibriumEnthalpyPhaseCurve {
    /// Admit a bounded, monotone solid-liquid enthalpy curve.
    ///
    /// At least two knots are required. Enthalpy must increase strictly;
    /// temperature and liquid fraction must not decrease. The first knot must
    /// be fully solid and the last fully liquid so that the curve owns both
    /// phase transitions rather than exposing a misleading partial taxonomy.
    pub fn try_new(
        material_card_identity: ContentHash,
        knots: Vec<EnthalpyPhaseKnot>,
    ) -> Result<Self, PhaseStateError> {
        validate_curve_common(material_card_identity, &knots)?;
        if knots[0].liquid_mass_fraction != 0.0
            || knots[knots.len() - 1].liquid_mass_fraction != 1.0
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "the curve must begin fully solid and end fully liquid",
            });
        }
        for pair in knots.windows(2) {
            if pair[1].temperature_k < pair[0].temperature_k {
                return Err(PhaseStateError::InvalidCurve {
                    what: "temperature must be nondecreasing with specific enthalpy",
                });
            }
            if pair[1].liquid_mass_fraction < pair[0].liquid_mass_fraction {
                return Err(PhaseStateError::InvalidCurve {
                    what: "liquid mass fraction must be nondecreasing with specific enthalpy",
                });
            }
        }
        Ok(Self::from_valid_knots(material_card_identity, knots))
    }

    /// Admit a bounded equilibrium enthalpy curve entirely in one explicit
    /// phase. Solid curves carry liquid fraction zero; liquid curves carry
    /// liquid fraction one. Temperature must rise strictly, so this ingress
    /// cannot hide a latent plateau or mixed transition.
    pub fn try_single_phase(
        material_card_identity: ContentHash,
        phase: SolidLiquidPhase,
        knots: Vec<EnthalpyPhaseKnot>,
    ) -> Result<Self, PhaseStateError> {
        validate_curve_common(material_card_identity, &knots)?;
        let liquid_mass_fraction = match phase {
            SolidLiquidPhase::Solid => 0.0,
            SolidLiquidPhase::Liquid => 1.0,
            SolidLiquidPhase::SolidLiquid => {
                return Err(PhaseStateError::InvalidCurve {
                    what: "a single-phase curve must declare Solid or Liquid",
                });
            }
        };
        if knots
            .iter()
            .any(|knot| knot.liquid_mass_fraction != liquid_mass_fraction)
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "every single-phase knot must carry the declared phase fraction",
            });
        }
        if knots
            .windows(2)
            .any(|pair| pair[1].temperature_k <= pair[0].temperature_k)
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "single-phase temperature knots must be strictly increasing",
            });
        }
        Ok(Self::from_valid_knots(material_card_identity, knots))
    }

    fn from_valid_knots(
        material_card_identity: ContentHash,
        knots: Vec<EnthalpyPhaseKnot>,
    ) -> Self {
        Self {
            material_card_identity,
            identity: phase_curve_identity(material_card_identity, &knots),
            knots,
        }
    }

    /// Material card that owns the supplied equilibrium data.
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card_identity
    }

    /// Complete phase-curve identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Admitted curve knots in ascending specific enthalpy.
    #[must_use]
    pub fn knots(&self) -> &[EnthalpyPhaseKnot] {
        &self.knots
    }

    /// Resolve equilibrium temperature, density, and phase fractions from
    /// specific enthalpy without extrapolation.
    pub fn state_at_specific_enthalpy(
        &self,
        specific_enthalpy_j_kg: f64,
    ) -> Result<EquilibriumPhaseState, PhaseStateError> {
        if !specific_enthalpy_j_kg.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnthalpy);
        }
        let lower = self.knots[0].specific_enthalpy_j_kg;
        let upper = self.knots[self.knots.len() - 1].specific_enthalpy_j_kg;
        if specific_enthalpy_j_kg < lower || specific_enthalpy_j_kg > upper {
            return Err(PhaseStateError::OutsideEnthalpyDomain {
                specific_enthalpy_j_kg,
                lower_j_kg: lower,
                upper_j_kg: upper,
            });
        }
        match self.knots.binary_search_by(|knot| {
            knot.specific_enthalpy_j_kg
                .total_cmp(&specific_enthalpy_j_kg)
        }) {
            Ok(index) => self.state_from_values(specific_enthalpy_j_kg, self.knots[index]),
            Err(upper_index) => {
                let lower_knot = self.knots[upper_index - 1];
                let upper_knot = self.knots[upper_index];
                let alpha = (specific_enthalpy_j_kg - lower_knot.specific_enthalpy_j_kg)
                    / (upper_knot.specific_enthalpy_j_kg - lower_knot.specific_enthalpy_j_kg);
                self.state_from_values(
                    specific_enthalpy_j_kg,
                    EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg,
                        temperature_k: lerp(
                            lower_knot.temperature_k,
                            upper_knot.temperature_k,
                            alpha,
                        ),
                        liquid_mass_fraction: lerp(
                            lower_knot.liquid_mass_fraction,
                            upper_knot.liquid_mass_fraction,
                            alpha,
                        ),
                        // Density is mass per volume, so interpolation on a
                        // mass-specific enthalpy coordinate is performed in
                        // specific volume and inverted. Linear density would
                        // violate additive mixture volume through a two-phase
                        // interval.
                        bulk_density_kg_m3: lerp(
                            lower_knot.bulk_density_kg_m3.recip(),
                            upper_knot.bulk_density_kg_m3.recip(),
                            alpha,
                        )
                        .recip(),
                    },
                )
            }
        }
    }

    /// Apply a signed specific-energy increment and resolve the resulting
    /// phase state. No heat source, boundary flux, or work term is invented;
    /// the caller owns that energy balance.
    pub fn advance_specific_energy(
        &self,
        current: EquilibriumPhaseState,
        net_specific_energy_j_kg: f64,
    ) -> Result<EquilibriumPhaseState, PhaseStateError> {
        if current.phase_curve_identity() != self.identity {
            return Err(PhaseStateError::StateCurveMismatch);
        }
        if !net_specific_energy_j_kg.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnergyIncrement);
        }
        let next = current.specific_enthalpy_j_kg() + net_specific_energy_j_kg;
        if !next.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnthalpy);
        }
        self.state_at_specific_enthalpy(next)
    }

    fn state_from_values(
        &self,
        specific_enthalpy_j_kg: f64,
        values: EnthalpyPhaseKnot,
    ) -> Result<EquilibriumPhaseState, PhaseStateError> {
        let liquid_mass_fraction = values.liquid_mass_fraction.clamp(0.0, 1.0);
        let solid_mass_fraction = 1.0 - liquid_mass_fraction;
        let phase = if liquid_mass_fraction == 0.0 {
            SolidLiquidPhase::Solid
        } else if liquid_mass_fraction == 1.0 {
            SolidLiquidPhase::Liquid
        } else {
            SolidLiquidPhase::SolidLiquid
        };
        let mut hasher = DomainHasher::new(EQUILIBRIUM_PHASE_STATE_IDENTITY_DOMAIN);
        hasher.update(self.identity.as_bytes());
        for value in [
            specific_enthalpy_j_kg,
            values.temperature_k,
            solid_mass_fraction,
            liquid_mass_fraction,
            values.bulk_density_kg_m3,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
        hasher.update(&[match phase {
            SolidLiquidPhase::Solid => 0,
            SolidLiquidPhase::SolidLiquid => 1,
            SolidLiquidPhase::Liquid => 2,
        }]);
        Ok(EquilibriumPhaseState {
            specific_enthalpy_j_kg,
            temperature_k: values.temperature_k,
            solid_mass_fraction,
            liquid_mass_fraction,
            bulk_density_kg_m3: values.bulk_density_kg_m3,
            phase,
            material_card_identity: self.material_card_identity,
            phase_curve_identity: self.identity,
            identity: hasher.finalize(),
        })
    }
}

fn validate_curve_common(
    material_card_identity: ContentHash,
    knots: &[EnthalpyPhaseKnot],
) -> Result<(), PhaseStateError> {
    if material_card_identity == ContentHash([0; 32]) {
        return Err(PhaseStateError::InvalidCurve {
            what: "material-card identity must not be zero",
        });
    }
    if knots.len() < 2 {
        return Err(PhaseStateError::InvalidCurve {
            what: "an equilibrium phase curve needs at least two knots",
        });
    }
    for knot in knots {
        if !(knot.specific_enthalpy_j_kg.is_finite()
            && knot.temperature_k > 0.0
            && knot.temperature_k.is_finite()
            && (0.0..=1.0).contains(&knot.liquid_mass_fraction)
            && knot.bulk_density_kg_m3 > 0.0
            && knot.bulk_density_kg_m3.is_finite())
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "every knot needs finite enthalpy, positive temperature/density, and liquid fraction in [0,1]",
            });
        }
    }
    for pair in knots.windows(2) {
        if pair[1].specific_enthalpy_j_kg <= pair[0].specific_enthalpy_j_kg
            || !(pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg).is_finite()
        {
            return Err(PhaseStateError::InvalidCurve {
                what: "specific enthalpy knots must be strictly increasing",
            });
        }
    }
    Ok(())
}

/// Typed refusal from equilibrium phase-state admission or evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum PhaseStateError {
    /// The source-provided curve violates a physical or ordering invariant.
    InvalidCurve {
        /// Failed curve invariant.
        what: &'static str,
    },
    /// A queried specific enthalpy was not finite.
    NonFiniteSpecificEnthalpy,
    /// A signed energy increment was not finite.
    NonFiniteSpecificEnergyIncrement,
    /// Evaluation would extrapolate beyond the source-provided curve.
    OutsideEnthalpyDomain {
        /// Refused query [J/kg].
        specific_enthalpy_j_kg: f64,
        /// Admitted lower endpoint [J/kg].
        lower_j_kg: f64,
        /// Admitted upper endpoint [J/kg].
        upper_j_kg: f64,
    },
    /// A state from a different phase curve was supplied to an update.
    StateCurveMismatch,
}

impl fmt::Display for PhaseStateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for PhaseStateError {}

fn phase_curve_identity(
    material_card_identity: ContentHash,
    knots: &[EnthalpyPhaseKnot],
) -> ContentHash {
    let mut hasher = DomainHasher::new(EQUILIBRIUM_ENTHALPY_PHASE_IDENTITY_DOMAIN);
    hasher.update(material_card_identity.as_bytes());
    hasher.update(&u64::try_from(knots.len()).unwrap_or(u64::MAX).to_le_bytes());
    for knot in knots {
        for value in [
            knot.specific_enthalpy_j_kg,
            knot.temperature_k,
            knot.liquid_mass_fraction,
            knot.bulk_density_kg_m3,
        ] {
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

fn lerp(start: f64, end: f64, alpha: f64) -> f64 {
    (end - start).mul_add(alpha, start)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plateau_curve() -> EquilibriumEnthalpyPhaseCurve {
        EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x51; 32]),
            vec![
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 0.0,
                    temperature_k: 300.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 11_300.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 30_000.0,
                    temperature_k: 600.0,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 11_100.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 55_000.0,
                    temperature_k: 600.0,
                    liquid_mass_fraction: 1.0,
                    bulk_density_kg_m3: 10_600.0,
                },
                EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: 95_000.0,
                    temperature_k: 800.0,
                    liquid_mass_fraction: 1.0,
                    bulk_density_kg_m3: 10_300.0,
                },
            ],
        )
        .unwrap()
    }

    #[test]
    fn g0_enthalpy_crosses_an_isothermal_phase_plateau_without_name_switches() {
        let curve = plateau_curve();
        let solid = curve.state_at_specific_enthalpy(20_000.0).unwrap();
        assert_eq!(solid.phase(), SolidLiquidPhase::Solid);
        assert_eq!(solid.temperature_k().to_bits(), 500.0_f64.to_bits());

        let half_melted = curve.advance_specific_energy(solid, 22_500.0).unwrap();
        assert_eq!(half_melted.phase(), SolidLiquidPhase::SolidLiquid);
        assert_eq!(half_melted.temperature_k().to_bits(), 600.0_f64.to_bits());
        assert_eq!(
            half_melted.solid_mass_fraction().to_bits(),
            0.5_f64.to_bits()
        );
        assert_eq!(
            half_melted.liquid_mass_fraction().to_bits(),
            0.5_f64.to_bits()
        );
        assert!(!half_melted.is_fully_solid());

        let liquid = curve
            .advance_specific_energy(half_melted, 32_500.0)
            .unwrap();
        assert_eq!(liquid.phase(), SolidLiquidPhase::Liquid);
        assert_eq!(liquid.temperature_k().to_bits(), 700.0_f64.to_bits());
        assert_eq!(liquid.liquid_mass_fraction().to_bits(), 1.0_f64.to_bits());
    }

    #[test]
    fn g0_phase_curve_refuses_extrapolation_and_foreign_state() {
        let curve = plateau_curve();
        assert!(matches!(
            curve.state_at_specific_enthalpy(95_000.0 + f64::EPSILON * 95_000.0),
            Err(PhaseStateError::OutsideEnthalpyDomain { .. })
        ));
        let other =
            EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x52; 32]), curve.knots().to_vec())
                .unwrap();
        let state = other.state_at_specific_enthalpy(40_000.0).unwrap();
        assert_eq!(
            curve.advance_specific_energy(state, 1.0),
            Err(PhaseStateError::StateCurveMismatch)
        );
    }

    #[test]
    fn g0_phase_curve_refuses_nonmonotone_or_partial_taxonomy() {
        let mut knots = plateau_curve().knots().to_vec();
        knots[3].liquid_mass_fraction = 0.25;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x53; 32]), knots),
            Err(PhaseStateError::InvalidCurve { .. })
        ));

        let mut knots = plateau_curve().knots().to_vec();
        knots[2].temperature_k = 599.0;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x54; 32]), knots),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
    }

    fn single_phase_knots(liquid_mass_fraction: f64) -> Vec<EnthalpyPhaseKnot> {
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 10.0,
                temperature_k: 300.0,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1_000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 110.0,
                temperature_k: 350.0,
                liquid_mass_fraction,
                bulk_density_kg_m3: 900.0,
            },
        ]
    }

    #[test]
    fn g0_single_phase_enthalpy_curves_interpolate_without_phase_transition() {
        let solid = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x61; 32]),
            SolidLiquidPhase::Solid,
            single_phase_knots(0.0),
        )
        .unwrap();
        let solid_midpoint = solid.state_at_specific_enthalpy(60.0).unwrap();
        assert_eq!(solid_midpoint.phase(), SolidLiquidPhase::Solid);
        assert_eq!(
            solid_midpoint.liquid_mass_fraction().to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            solid_midpoint.solid_mass_fraction().to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(
            solid_midpoint.temperature_k().to_bits(),
            325.0_f64.to_bits()
        );
        assert!((solid_midpoint.bulk_density_kg_m3() - 947.368_421_052_631_6).abs() < 1e-12);

        let liquid = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x62; 32]),
            SolidLiquidPhase::Liquid,
            single_phase_knots(1.0),
        )
        .unwrap();
        let liquid_midpoint = liquid.state_at_specific_enthalpy(60.0).unwrap();
        assert_eq!(liquid_midpoint.phase(), SolidLiquidPhase::Liquid);
        assert_eq!(
            liquid_midpoint.liquid_mass_fraction().to_bits(),
            1.0_f64.to_bits()
        );
        assert_eq!(
            liquid_midpoint.solid_mass_fraction().to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            liquid_midpoint.temperature_k().to_bits(),
            325.0_f64.to_bits()
        );
    }

    #[test]
    fn g0_single_phase_enthalpy_curves_refuse_transition_plateau_and_extrapolation() {
        let identity = ContentHash([0x63; 32]);
        let mut mixed_fraction = single_phase_knots(0.0);
        mixed_fraction[1].liquid_mass_fraction = 0.5;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                identity,
                SolidLiquidPhase::Solid,
                mixed_fraction,
            ),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                identity,
                SolidLiquidPhase::SolidLiquid,
                single_phase_knots(0.0),
            ),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
        let mut plateau = single_phase_knots(0.0);
        plateau[1].temperature_k = plateau[0].temperature_k;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                identity,
                SolidLiquidPhase::Solid,
                plateau
            ),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
        let mut cooling = single_phase_knots(1.0);
        cooling[1].temperature_k = 299.0;
        assert!(matches!(
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                identity,
                SolidLiquidPhase::Liquid,
                cooling
            ),
            Err(PhaseStateError::InvalidCurve { .. })
        ));
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            identity,
            SolidLiquidPhase::Solid,
            single_phase_knots(0.0),
        )
        .unwrap();
        assert!(matches!(
            curve.state_at_specific_enthalpy(9.0),
            Err(PhaseStateError::OutsideEnthalpyDomain { .. })
        ));
    }
}
