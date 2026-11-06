//! Source-backed entropy-constrained waves in a finitely strained solid.
//!
//! The native finite thermal jet differentiates the original full compressible
//! free energy at fixed TOTAL F. The same-card caloric chart supplies its own
//! local capacity; original reference density converts that specific capacity.
//! This local reversible response reuses the existing acoustic tensor, signed
//! modes and energy velocities. It supplies no heat evolution or attenuation.
//!
//! Expansion paths and caloric charts need a smooth two-sided neighborhood.
//! Source kinks, phase changes and mean-coefficient point evidence cannot
//! silently supply missing derivatives. The original stress-free zero-pressure
//! caloric convention is explicit; arbitrary enthalpy is not converted to u.
//! Exact-incompressible and retained-volume mixed bodies require their own
//! constrained acoustic linearization and are outside this per-material view.
//! Native constitutive activation and switching tangent conventions are retained;
//! smooth thermal/caloric paths do not certify two-sided mechanical derivatives
//! at a fiber-activation boundary.

mod coefficients;

use super::conservative::ConservativeThermalField;
use super::{
    ExpansionThermalFrame, ExpansionThermalReference, FiniteSolidMaterial, ThermalReference,
};
pub use crate::solid_enthalpy::finite_plastic::CaloricReference;
use fs_material::MaterialError;
use fs_material::hyper::acoustic::DirectionalHyperelasticWaves;
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::thermoelastic::source::{
    NeoHookeanThermalSourceError, ResolvedNeoHookeanThermalCoefficients,
};
use fs_material::hyper::thermoelastic::{
    AdiabaticHyperelasticState, FiniteThermoelasticState, ThermalReferenceJet,
};
use fs_material::phase::{EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve};
use fs_material::thermal::transport::ThermalTransportError;
use fs_qty::semantic::QuantityKind;

#[derive(Debug)]
pub enum FiniteThermoelasticSourceError {
    Invalid(&'static str),
    Material(MaterialError),
    Caloric(PhaseStateError),
    Expansion(ThermalTransportError),
    Coefficients(NeoHookeanThermalSourceError),
}

impl core::fmt::Display for FiniteThermoelasticSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite thermoelastic source refused: {self:?}")
    }
}
impl core::error::Error for FiniteThermoelasticSourceError {}

fn invalid(message: &'static str) -> FiniteThermoelasticSourceError {
    FiniteThermoelasticSourceError::Invalid(message)
}

// These comparisons admit arithmetic roundoff at a redundant source knot;
// they do not smooth a measured slope change or create a temperature interval.
fn same_derivative(a: f64, b: f64) -> bool {
    a == b || (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

fn expansion_jet(
    source: &ExpansionThermalReference,
    before_thermal: ReferenceDeformation,
) -> Result<ThermalReferenceJet, FiniteThermoelasticSourceError> {
    let path = source.source();
    let temperature = source.geometry().temperature_k();
    let [lower, upper] = path.temperature_span_k();
    if temperature <= lower || temperature >= upper {
        return Err(invalid(
            "adiabatic waves require a two-sided expansion-source temperature neighborhood",
        ));
    }
    let tangent = path
        .at_temperature_with_tangent(temperature)
        .map_err(FiniteThermoelasticSourceError::Expansion)?;
    let slopes = tangent.principal_stretch_derivatives_per_k();
    let mut previous = None;
    for knot in path.temperatures_k() {
        if knot == temperature {
            // At the preceding knot the native right-interval derivative is
            // exactly the left derivative at this knot. No differencing or
            // resampling of mechanical energy is used to infer a source slope.
            let left = path
                .at_temperature_with_tangent(
                    previous.ok_or_else(|| {
                        invalid("a two-sided expansion-source interval is required")
                    })?,
                )
                .map_err(FiniteThermoelasticSourceError::Expansion)?
                .principal_stretch_derivatives_per_k();
            if left
                .into_iter()
                .zip(slopes)
                .any(|(a, b)| !same_derivative(a, b))
            {
                return Err(invalid(
                    "a slope-changing expansion knot has no unique thermoelastic linearization",
                ));
            }
            break;
        }
        if knot > temperature {
            break;
        }
        previous = Some(knot);
    }
    let q = match source.frame() {
        ExpansionThermalFrame::Isotropic => [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ExpansionThermalFrame::FixedPrincipal {
            principal_to_mesh, ..
        } => principal_to_mesh,
    };
    let thermal_derivative: [f64; 9] = core::array::from_fn(|ij| {
        (0..3)
            .map(|k| q[ij / 3][k] * slopes[k] * q[ij % 3][k])
            .sum()
    });
    let before = before_thermal.deformation_gradient();
    let derivative = core::array::from_fn(|ij| {
        (0..3)
            .map(|k| thermal_derivative[3 * (ij / 3) + k] * before[3 * k + ij % 3])
            .sum()
    });
    let reference = before_thermal
        .followed_by(*source.reference_deformation())
        .map_err(FiniteThermoelasticSourceError::Material)?;
    // The original engineering-length path is piecewise linear with fixed
    // principal axes. At the admitted smooth point G_TT is exactly zero;
    // multiplying by a held manufacturing map preserves that derivative.
    ThermalReferenceJet::try_new(reference, derivative, [0.0; 9])
        .map_err(FiniteThermoelasticSourceError::Material)
}

fn caloric_capacity(
    material: &FiniteSolidMaterial,
    caloric: &ResolvedHeatCapacityCurve,
    h: f64,
    temperature: f64,
    reference: CaloricReference,
) -> Result<(EquilibriumPhaseState, f64), FiniteThermoelasticSourceError> {
    match reference {
        CaloricReference::StressFreeZeroPressure => {}
    }
    let curve = caloric.curve();
    let phase = curve
        .state_at_specific_enthalpy(h)
        .map_err(FiniteThermoelasticSourceError::Caloric)?;
    let context = material.law().context();
    if curve.material_card_identity() != material.density().card_identity()
        || !phase.is_fully_solid()
        || phase.temperature_k() != temperature
    {
        return Err(invalid(
            "adiabatic elasticity requires the same solid material card and actual caloric temperature",
        ));
    }
    if context.axis_quantities().iter().any(|(axis, quantity)| {
        quantity
            .semantic_type()
            .is_some_and(|kind| kind.kind() == QuantityKind::Pressure)
            && context.axes().get(axis).is_some_and(|value| *value != 0.0)
    }) {
        return Err(invalid(
            "the zero-pressure caloric reference cannot use a nonzero pressure-conditioned context",
        ));
    }
    if material
        .law()
        .model_card()
        .validity
        .bound("T")
        .is_some_and(|(low, high)| temperature <= low || temperature >= high)
    {
        return Err(invalid(
            "adiabatic waves require a two-sided constitutive source temperature neighborhood",
        ));
    }
    let receipts = caloric.receipts();
    let density_claim = material.density().properties()[0].answer().receipt.selected;
    if receipts.len() < 4
        || receipts.len() % 2 != 0
        || receipts.chunks_exact(2).any(|pair| {
            pair[0].selected != receipts[0].selected || pair[1].selected != density_claim
        })
        || receipts.iter().any(|receipt| {
            receipt.axis_quantities != *context.axis_quantities()
                || !receipt
                    .query_point
                    .iter()
                    .filter(|(axis, _)| axis != "T")
                    .map(|(axis, value)| (axis, value))
                    .eq(context
                        .axes()
                        .iter()
                        .filter(|(axis, _)| axis.as_str() != "T"))
        })
    {
        return Err(invalid(
            "caloric data must retain the original density claim, pinned capacity source and fixed typed context",
        ));
    }
    let knots = curve.knots();
    if knots.len() < 2
        || h <= knots[0].specific_enthalpy_j_kg
        || h >= knots[knots.len() - 1].specific_enthalpy_j_kg
        || knots.iter().any(|knot| knot.liquid_mass_fraction != 0.0)
    {
        return Err(invalid(
            "adiabatic solid waves require a two-sided fully solid caloric neighborhood",
        ));
    }
    let upper = knots.partition_point(|knot| knot.specific_enthalpy_j_kg <= h);
    let first = if knots[upper - 1].specific_enthalpy_j_kg == h {
        upper - 2
    } else {
        upper - 1
    };
    let mut previous = None;
    for (offset, pair) in knots[first..=upper].windows(2).enumerate() {
        let capacity = curve
            .sensible_heat_capacities_on_interval(first + offset)
            .map_or_else(
                || {
                    (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                        / (pair[1].temperature_k - pair[0].temperature_k)
                },
                |[lower, upper]| {
                    let fraction = (temperature - pair[0].temperature_k)
                        / (pair[1].temperature_k - pair[0].temperature_k);
                    lower * (1.0 - fraction) + upper * fraction
                },
            );
        if !capacity.is_finite()
            || capacity <= 0.0
            || previous.is_some_and(|left| !same_derivative(left, capacity))
        {
            return Err(invalid(
                "a nonpositive capacity or slope-changing caloric knot has no unique adiabatic linearization",
            ));
        }
        previous = Some(capacity);
    }
    let derivative = curve
        .temperature_derivative_at_specific_enthalpy(h)
        .map_err(FiniteThermoelasticSourceError::Caloric)?;
    let capacity = material.reference_density_kg_m3() / derivative;
    if !capacity.is_finite() || capacity <= 0.0 {
        return Err(invalid(
            "positive finite caloric capacity per original volume required",
        ));
    }
    Ok((phase, capacity))
}

/// Same-source local thermodynamics and finite-strain isentropic wave response.
/// Caloric density metadata is retained as evidence; inertial density remains
/// this material's original reference mass divided by its original volume.
#[derive(Debug)]
pub struct FiniteSolidAdiabaticState<'a> {
    material: &'a FiniteSolidMaterial,
    caloric: &'a ResolvedHeatCapacityCurve,
    phase: EquilibriumPhaseState,
    before_thermal: ReferenceDeformation,
    native: AdiabaticHyperelasticState,
    coefficients: Option<&'a ResolvedNeoHookeanThermalCoefficients>,
}

impl FiniteSolidMaterial {
    /// Prepare a full compressible local response using G=G_thermal*G_before.
    /// The source's original static model parameters stay fixed; temperature
    /// dependence is supplied by its admitted engineering-length path.
    ///
    /// The caloric chart must declare a stress-free zero-pressure reference,
    /// match the source's actual temperature and retain its original density
    /// claim/context. The native law adds mechanical capacity exactly once and
    /// refuses a nonpositive total capacity. No full-interval capacity bound,
    /// dynamic heat closure or measured ultrasonic validity is implied.
    ///
    /// This uses the per-material reference density. A body's separately bound
    /// density field requires its own local inertial/caloric inventory.
    #[allow(clippy::too_many_arguments)]
    pub fn adiabatic_acoustic_state<'a>(
        &'a self,
        total_f: &[f64; 9],
        before_thermal: ReferenceDeformation,
        caloric: &'a ResolvedHeatCapacityCurve,
        specific_h_j_kg: f64,
        caloric_reference: CaloricReference,
    ) -> Result<FiniteSolidAdiabaticState<'a>, FiniteThermoelasticSourceError> {
        self.adiabatic_acoustic_state_inner(
            total_f,
            before_thermal,
            caloric,
            specific_h_j_kg,
            caloric_reference,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn adiabatic_acoustic_state_inner<'a>(
        &'a self,
        total_f: &[f64; 9],
        before_thermal: ReferenceDeformation,
        caloric: &'a ResolvedHeatCapacityCurve,
        specific_h_j_kg: f64,
        caloric_reference: CaloricReference,
        coefficients: Option<&'a ResolvedNeoHookeanThermalCoefficients>,
    ) -> Result<FiniteSolidAdiabaticState<'a>, FiniteThermoelasticSourceError> {
        let Some(ThermalReference::Expansion(expansion)) = self.thermal_reference() else {
            return Err(invalid(
                "finite thermoelastic source derivatives require a prepared expansion path; mean-coefficient endpoint evidence is insufficient",
            ));
        };
        let temperature = expansion.geometry().temperature_k();
        if coefficients.is_some_and(|source| {
            source.model() != self.law()
                || source.material().card_identity() != self.density().card_identity()
                || source.temperature_k() != temperature
        }) {
            return Err(invalid(
                "thermal coefficient derivatives must retain the exact original material/model/context at its actual temperature",
            ));
        }
        let jet = expansion_jet(expansion, before_thermal)?;
        let (phase, capacity) = caloric_capacity(
            self,
            caloric,
            specific_h_j_kg,
            temperature,
            caloric_reference,
        )?;
        let potential = match coefficients {
            Some(source) => source
                .thermoelastic_state(total_f, jet)
                .map_err(FiniteThermoelasticSourceError::Coefficients)?,
            None => FiniteThermoelasticState::try_new(self.law().law(), total_f, temperature, jet)
                .map_err(FiniteThermoelasticSourceError::Material)?,
        };
        let native = potential
            .adiabatic_acoustics(self.reference_density_kg_m3(), capacity)
            .map_err(FiniteThermoelasticSourceError::Material)?;
        Ok(FiniteSolidAdiabaticState {
            material: self,
            caloric,
            phase,
            before_thermal,
            native,
            coefficients,
        })
    }
}

impl ConservativeThermalField {
    /// Full compressible per-material response with this field's exact source
    /// and manufacturing map. This does not infer an accepted mixed body's
    /// retained volume, interval pressure or constrained acoustic projection.
    pub fn adiabatic_acoustic_state<'a>(
        &'a self,
        element: usize,
        total_f: &[f64; 9],
        caloric: &'a ResolvedHeatCapacityCurve,
        specific_h_j_kg: f64,
        caloric_reference: CaloricReference,
    ) -> Result<FiniteSolidAdiabaticState<'a>, FiniteThermoelasticSourceError> {
        let material = self
            .material_at(element)
            .ok_or_else(|| invalid("thermoelastic source element is outside the original field"))?;
        let before = self
            .reference_before_thermal_at(element)
            .ok_or_else(|| invalid("the original manufacturing map is absent"))?;
        material.adiabatic_acoustic_state(
            total_f,
            *before,
            caloric,
            specific_h_j_kg,
            caloric_reference,
        )
    }
}

impl<'a> FiniteSolidAdiabaticState<'a> {
    /// Original E(T)/nu(T) evidence when the coefficient-dependent path was
    /// explicitly selected; absent on the existing fixed-coefficient path.
    pub const fn thermal_coefficients(&self) -> Option<&'a ResolvedNeoHookeanThermalCoefficients> {
        self.coefficients
    }
    pub const fn material(&self) -> &'a FiniteSolidMaterial {
        self.material
    }
    pub const fn caloric_source(&self) -> &'a ResolvedHeatCapacityCurve {
        self.caloric
    }
    pub const fn caloric_state(&self) -> EquilibriumPhaseState {
        self.phase
    }
    pub const fn reference_before_thermal(&self) -> &ReferenceDeformation {
        &self.before_thermal
    }
    pub const fn native(&self) -> &AdiabaticHyperelasticState {
        &self.native
    }
    pub fn thermodynamics(&self) -> &FiniteThermoelasticState {
        self.native.source()
    }
    pub fn directional_waves(
        &self,
        direction: [f64; 3],
    ) -> Result<DirectionalHyperelasticWaves, MaterialError> {
        self.native.acoustic().directional_waves(direction)
    }
}
