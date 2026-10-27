//! Source-backed modified Arrhenius kinetics, with explicit normalization.
//!
//! `k(T) = A (T/T0)^b exp(-Ea/(R T))`: A is in mol/(m3 s), b is
//! dimensionless, signed effective Ea is in J/mol, and T0 is in kelvin.
//! Pressure activities give `r_f = k(T) product_i (p_i/p0)^alpha_i`.
//! Normalized concentrations give `r_f = k(T) product_i (c_i/c0)^alpha_i`.
//! Since `c_i = p_i/(R T)`, the original activity-based owner receives
//! `k_activity = k(T) (p0/(R T c0))^sum(alpha_i)`. Forward orders come from
//! the declared stoichiometry; the original reverse closure preserves Kp.
//!
//! Importing a concentration coefficient A_c requires the explicit conversion
//! `A = A_c c0^sum(alpha_i)` to mol, m, s and the selected T/T0 convention.
//! No third-body, falloff, surface kinetics or unlabeled database units are
//! inferred. Reference: <https://www.cantera.org/stable/reference/kinetics/rate-constants.html>.
//!
//! Coefficient claim IDs are pinned before evolution. Every explicit/implicit
//! evaluation uses the original gas query point and material source resolver.
//! The source bundle retains coefficients and the kinetic receipt retains the
//! effective scale. Domain/range refusals occur before coupled publication.
//! No new integrator, inventory, thermal source or gas constant is introduced.
use super::{
    ACTIVITY_RATE_SCALE_DIMS, DeclaredStoichiometricMassActionV1, GasReactionError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement, Temperature,
    invalid, positive,
};
use fs_matdb::MaterialCard;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, Pressure};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Prefactor in mol/(m3 s) for explicitly normalized activity/concentration.
pub const PREFACTOR_PROPERTY: &str = "reaction_arrhenius_prefactor";
/// Dimensionless exponent on T/T0.
pub const TEMPERATURE_EXPONENT_PROPERTY: &str = "reaction_arrhenius_temperature_exponent";
/// Signed effective activation energy in J/mol, not J/kmol or kelvin.
pub const ACTIVATION_ENERGY_PROPERTY: &str = "reaction_arrhenius_activation_energy";
/// Absolute temperature normalization T0 in kelvin.
pub const REFERENCE_TEMPERATURE_PROPERTY: &str = "reaction_arrhenius_reference_temperature";
/// Explicit concentration normalization c0 in mol/m3.
pub const REFERENCE_CONCENTRATION_PROPERTY: &str = "reaction_arrhenius_reference_concentration";
/// Molar energy dimensions J/mol.
pub const MOLAR_ACTIVATION_ENERGY_DIMS: Dims = Dims([2, 1, -2, 0, 0, -1]);
/// Molar concentration dimensions mol/m3.
pub const MOLAR_CONCENTRATION_DIMS: Dims = Dims([-3, 0, 0, 0, 0, 1]);

/// Explicit data-author convention, not an integration option.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ArrheniusRateBasis {
    /// The existing dimensionless ideal-gas activities p_i/p0.
    PressureActivity,
    /// Dimensionless molar concentrations c_i/c0, with material-supplied c0.
    NormalizedMolarConcentration,
}

pub(super) fn pin_coefficients(
    material: &MaterialCard,
    selection: &MaterialPropertySelection,
    basis: ArrheniusRateBasis,
) -> Result<MaterialPropertySelection, GasReactionError> {
    pin_requirements(material, selection, &requirements(basis)?)
}

// One pinning policy for all kinetic source properties. Collider laws may
// extend the required set, but cannot select fresh claims during evolution.
pub(super) fn pin_requirements(
    material: &MaterialCard,
    selection: &MaterialPropertySelection,
    required: &[ScalarPropertyRequirement],
) -> Result<MaterialPropertySelection, GasReactionError> {
    match selection {
        // Full keys, coverage, source membership and point validity remain the
        // original shared resolver's responsibility on every evaluation.
        MaterialPropertySelection::PinnedByProperty(_) => Ok(selection.clone()),
        MaterialPropertySelection::PreferObservationBacked => Err(invalid(
            "kinetic source properties require explicit claim pins after observation selection",
        )),
        MaterialPropertySelection::SingleClaimOnly => {
            let mut pins = Vec::with_capacity(required.len());
            for requirement in required {
                let mut claims = material
                    .claims()
                    .claims_ordered()
                    .filter(|(_, claim)| claim.key.name() == requirement.name());
                let Some((id, _)) = claims.next() else {
                    return Err(invalid(
                        "kinetic material is missing a required property claim",
                    ));
                };
                if claims.next().is_some() {
                    return Err(invalid(
                        "kinetic property has multiple claims; explicitly pin one before evolution",
                    ));
                }
                pins.push((requirement.name().to_owned(), id));
            }
            Ok(MaterialPropertySelection::PinnedByProperty(pins))
        }
    }
}

pub(super) fn requirements(
    basis: ArrheniusRateBasis,
) -> Result<Vec<ScalarPropertyRequirement>, GasReactionError> {
    use ScalarAdmissibility::{Finite, StrictlyPositive};
    let mut required = Vec::with_capacity(5);
    for (name, dims, domain) in [
        (
            PREFACTOR_PROPERTY,
            ACTIVITY_RATE_SCALE_DIMS,
            StrictlyPositive,
        ),
        (TEMPERATURE_EXPONENT_PROPERTY, Dims([0; 6]), Finite),
        (
            ACTIVATION_ENERGY_PROPERTY,
            MOLAR_ACTIVATION_ENERGY_DIMS,
            Finite,
        ),
        (
            REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            StrictlyPositive,
        ),
    ] {
        required.push(ScalarPropertyRequirement::try_new(name, dims, domain)?);
    }
    if basis == ArrheniusRateBasis::NormalizedMolarConcentration {
        required.push(ScalarPropertyRequirement::try_new(
            REFERENCE_CONCENTRATION_PROPERTY,
            MOLAR_CONCENTRATION_DIMS,
            StrictlyPositive,
        )?);
    }
    Ok(required)
}

pub(super) fn activity_scale(
    resolved: &ResolvedMaterialStatePoint,
    temperature: Temperature,
    reference_pressure: Pressure,
    basis: ArrheniusRateBasis,
    kinetics: &DeclaredStoichiometricMassActionV1,
) -> Result<f64, GasReactionError> {
    let property = |name| {
        resolved
            .property(name)
            .map(|p| p.value_si())
            .ok_or_else(|| invalid("resolved material omitted a required Arrhenius coefficient"))
    };
    let t = temperature.value();
    let a = property(PREFACTOR_PROPERTY)?;
    let b = property(TEMPERATURE_EXPONENT_PROPERTY)?;
    let ea = property(ACTIVATION_ENERGY_PROPERTY)?;
    let t0 = property(REFERENCE_TEMPERATURE_PROPERTY)?;
    if !positive(t) || !positive(reference_pressure.value()) {
        return Err(invalid(
            "Arrhenius evaluation needs positive absolute temperature and reference pressure",
        ));
    }
    // Log ratios avoid avoidable intermediate overflow/underflow. Neither
    // signed effective Ea nor b is clipped; a positive scale must survive.
    let mut log_scale =
        fs_math::det::ln(a) + b * (fs_math::det::ln(t) - fs_math::det::ln(t0)) - (ea / R) / t;
    if basis == ArrheniusRateBasis::NormalizedMolarConcentration {
        let c0 = property(REFERENCE_CONCENTRATION_PROPERTY)?;
        let order: f64 = kinetics
            .reverse_rate_closure()
            .equilibrium()
            .terms()
            .iter()
            .filter(|term| term.coefficient() < 0)
            .map(|term| -(term.coefficient() as f64))
            .sum();
        if !positive(order) {
            return Err(invalid(
                "Arrhenius concentration conversion needs represented forward kinetic order",
            ));
        }
        log_scale += order
            * (fs_math::det::ln(reference_pressure.value())
                - fs_math::det::ln(R)
                - fs_math::det::ln(t)
                - fs_math::det::ln(c0));
    }
    if !log_scale.is_finite() {
        return Err(invalid("Arrhenius log scale is unrepresentable"));
    }
    let scale = fs_math::det::exp(log_scale);
    if !positive(scale) {
        return Err(invalid(
            "Arrhenius scale underflows or overflows the kinetic quantity domain",
        ));
    }
    Ok(scale)
}
