//! Source-bound thermoelastic storage for the existing live uniform string.
//!
//! The fixed support length, section and reference mass retain the native sine
//! basis. The original expansion path supplies eigenstrain; a scalar or
//! continuous piecewise-affine alpha and a globally affine Young's modulus
//! supply the temperature derivatives of axial and transverse Helmholtz
//! energies. Alpha may change slope at a source knot: entropy stays continuous,
//! while the one-sided heat capacities may differ. No fitted
//! frequency correction or second acoustic/thermal solver is introduced.
//!
//! Axial energy uses A*L_ref, so its derivative with respect to loaded length
//! is the existing EA*elastic_strain tension. Caloric mass and transverse modal
//! normalization retain the existing loaded geometry. This is the declared
//! small-strain, fixed-section model, without transverse contraction or volume
//! work. Internal energy is Psi+T*S; the kinetic energy is temperature independent.

use super::{LinearMaterialStringRuntime, ThermalMaterialStringRuntime, ThermalStringError};
use crate::string_specimen::{ResolvedStringSpecimen, StringPrestress,
    with_uniform_circular_thermal_extension};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState,
    MAX_HEAT_CAPACITY_ENTHALPY_KNOTS};
use fs_material::state_point::{IntegratedIsotropicThermalExpansion,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    affine_temperature_slope, integrate_isotropic_thermal_expansion,
    resolve_material_state_point};
use fs_matdb::{MaterialCard, PropertyValue};
use std::sync::Arc;

/// Bounded reversible-heat solve on the original caloric/transport owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoelasticStringSolve {
    /// Caller-supplied signed bracket for reversible heat into the caloric
    /// store [J], containing zero and the actual coupled residual's root.
    pub reversible_heat_bracket_j: [f64; 2],
    /// Absolute coupled heat/energy residual allowance [J]. This is a numerical
    /// solve tolerance, not a material or temporal-discretization error bound.
    pub energy_tolerance_j: f64,
    /// Hard maximum number of scalar iterations in 1..=128; lack of convergence refuses.
    pub maximum_iterations: usize,
}

/// The sourced law and loading datum, independent of changing modal state.
#[derive(Clone)]
pub(crate) struct ThermoelasticString {
    controls: ThermoelasticStringSolve,
    expansion: IntegratedIsotropicThermalExpansion,
    reference_length_m: f64,
    strain_limit: f64,
    initial_temperature_k: f64,
    initial_modulus_pa: f64,
    initial_elastic_strain: f64,
    modulus_slope_pa_k: f64,
    expansion_per_k: f64,
    // None retains the exact constant-alpha capacity admission. A varying
    // original table stores only conservative capacity bounds, never another
    // expansion law or a resampled eigenstrain.
    expansion_capacity: Option<Arc<[ExpansionCapacityBound]>>,
}

/// Prepared upper bound on Psi_TT over one original CTE/caloric intersection.
/// The modal gradient is the only changing factor; the reference inventory,
/// source slopes and admitted small-strain envelope are immutable.
#[derive(Clone, Copy)]
struct ExpansionCapacityBound {
    upper_temperature_k: f64,
    caloric_capacity_j_k: f64,
    axial_curvature_j_k2: f64,
    modal_curvature_per_gradient: f64,
}

/// Free-energy and entropy terms used by the shared coupled first-law balance.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ThermoelasticEnergy {
    pub axial_free_energy_j: f64,
    pub modal_free_energy_j: f64,
    /// Axial and transverse contributions, excluding caloric entropy [J/K].
    pub entropy_j_k: f64,
}

impl ThermalMaterialStringRuntime {
    /// Attach fixed-support thermoelastic feedback to the existing live string.
    ///
    /// The specimen must retain its original sourced thermal-expansion path.
    /// Scalar or continuous piecewise-affine expansion and globally affine
    /// Young's modulus are admitted on the same immutable card. A varying
    /// expansion source must cover the original solid caloric interval; a
    /// conservative whole-interval caloric-plus-elastic capacity guard includes
    /// its original slopes and the actual retained modal amplitudes.
    /// Density and the section remain fixed. Heating
    /// changes tension, modes and diagnostic pressure without resetting motion.
    /// The shared thermal, electrical and finite-gas paths account for axial
    /// storage and the entropy of both axial and transverse elastic energies.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_thermoelastic(
        cx: &fs_exec::Cx<'_>,
        mechanical: LinearMaterialStringRuntime,
        card: MaterialCard,
        curve: EquilibriumEnthalpyPhaseCurve,
        specific_enthalpy_j_kg: f64,
        maximum_temperature_increment_k: f64,
        controls: ThermoelasticStringSolve,
    ) -> Result<Self, ThermalStringError> {
        let law = ThermoelasticString::try_new(cx, &mechanical, &card, &curve, controls)?;
        Self::try_new_impl(cx, mechanical, card, curve, specific_enthalpy_j_kg,
            maximum_temperature_increment_k, Some(law))
    }
}

impl ThermoelasticString {
    pub(crate) fn try_new(
        cx: &fs_exec::Cx<'_>,
        mechanical: &LinearMaterialStringRuntime,
        card: &MaterialCard,
        curve: &EquilibriumEnthalpyPhaseCurve,
        controls: ThermoelasticStringSolve,
    ) -> Result<Self, ThermalStringError> {
        check(cx)?;
        let [lower, upper] = controls.reversible_heat_bracket_j;
        if !lower.is_finite() || !upper.is_finite() || lower >= upper || lower > 0.0 || upper < 0.0
            || !positive(controls.energy_tolerance_j) || controls.maximum_iterations == 0
            || controls.maximum_iterations > 128
        {
            return Err(admit("thermoelastic heat solve needs a finite ordered bracket containing zero, positive tolerance and 1..=128 iterations"));
        }
        let specimen = mechanical.specimen();
        let StringPrestress::FixedThermalExtension {
            reference_stress_free_length_m, free_thermal_strain, linear_strain_limit,
        } = specimen.prestress() else {
            return Err(admit("thermoelastic feedback requires fixed-support thermal extension"));
        };
        let expansion = specimen.thermal_expansion().ok_or_else(||
            admit("thermoelastic feedback requires the original sourced expansion path"))?;
        let material = specimen.material();
        let initial_temperature_k = material.query().axes().get("T").copied()
            .ok_or_else(|| admit("thermoelastic source requires absolute temperature"))?;
        if material.card_identity() != card.content_hash()
            || expansion.current().resolved().card_identity() != material.card_identity()
            || expansion.current().resolved().query_point() != material.query_point()
            || expansion.current().resolved().axis_quantities() != material.axis_quantities()
            || expansion.free_linear_strain() != free_thermal_strain
        {
            return Err(admit("elasticity and expansion must retain the same card and physical state"));
        }
        let alpha = expansion.current().resolved()
            .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .ok_or_else(|| admit("retained expansion source lacks its coefficient"))?;
        let modulus = material.property("young_modulus")
            .ok_or_else(|| admit("thermoelastic string needs Young's modulus"))?;
        let claim = card.claims().claim(modulus.answer().receipt.selected)
            .ok_or_else(|| admit("selected elastic claim is absent"))?;
        let modulus_slope_pa_k = match &claim.value {
            PropertyValue::Scalar { .. } => 0.0,
            PropertyValue::Curve { abscissa, knots, .. }
                if abscissa == "T" && claim.interpolation == fs_matdb::InterpolationPolicy::LinearInside
                    && knots.len() <= MAX_HEAT_CAPACITY_ENTHALPY_KNOTS =>
            {
                affine_temperature_slope(knots, initial_temperature_k, modulus.value_si())
                    .map_err(ThermalStringError::Admission)?
            }
            _ => return Err(admit("Young's modulus must be scalar or one bounded globally affine temperature curve")),
        };
        let initial_elastic_strain = finite(
            (specimen.string().length_m-reference_stress_free_length_m)
                / reference_stress_free_length_m-free_thermal_strain)?;
        let mut law = Self { controls, expansion: expansion.clone(),
            reference_length_m: reference_stress_free_length_m, strain_limit: linear_strain_limit,
            initial_temperature_k, initial_modulus_pa: modulus.value_si(),
            initial_elastic_strain, modulus_slope_pa_k, expansion_per_k: alpha.value_si(),
            expansion_capacity: None };
        law.expansion_capacity = law.prepare_expansion_capacity(cx, mechanical, card, curve)?;
        law.admit_capacity(cx, mechanical, curve)?;
        law.energy(mechanical, initial_temperature_k)?;
        check(cx)?;
        Ok(law)
    }

    pub(crate) const fn controls(&self) -> ThermoelasticStringSolve { self.controls }

    /// Resolve through the existing material, expansion and specimen owners.
    pub(crate) fn specimen_at(
        &self,
        runtime: &ThermalMaterialStringRuntime,
        thermal: EquilibriumPhaseState,
    ) -> Result<ResolvedStringSpecimen, ThermalStringError> {
        let point = runtime.point.clone().with("T", thermal.temperature_k())
            .map_err(ThermalStringError::Query)?;
        let state = resolve_material_state_point(&runtime.card, &point,
            &runtime.requirements, runtime.selection.clone()).map_err(ThermalStringError::Material)?;
        let old = runtime.mechanical.specimen();
        // Source support is a path property. Valid endpoints cannot bridge a gap.
        let property = old.material().property("young_modulus")
            .ok_or_else(|| admit("thermoelastic string needs Young's modulus"))?;
        let selected = property.answer().receipt.selected;
        let claim = runtime.card.claims().claim(selected)
            .ok_or_else(|| admit("selected elastic claim is absent"))?;
        let from = old.material().query();
        let (lower, upper) = if from.axes()["T"] <= point.axes()["T"] {
            (&from, &point)
        } else { (&point, &from) };
        runtime.card.claims().query_line_typed(&claim.key, "T", lower, upper,
            fs_matdb::ClaimSelection::Pinned(selected))
            .map_err(|_| admit("elastic source does not support the complete thermal path"))?;
        let expansion = integrate_isotropic_thermal_expansion(&runtime.card,
            &self.expansion.reference().resolved().query(), &point,
            MaterialPropertySelection::PinnedByProperty(vec![(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY.to_owned(),
                self.expansion.selected_claim(),
            )])).map_err(ThermalStringError::Material)?;
        with_uniform_circular_thermal_extension(old.string(), old.geometry_constraint(),
            &state, &expansion, self.reference_length_m, self.strain_limit)
            .and_then(|specimen| specimen.with_kelvin_voigt_bending_loss())
            .map_err(ThermalStringError::Acoustic)
    }

    /// Evaluate Psi and -dPsi/dT at fixed modal coordinates and support length.
    pub(crate) fn energy(
        &self,
        mechanical: &LinearMaterialStringRuntime,
        temperature_k: f64,
    ) -> Result<ThermoelasticEnergy, ThermalStringError> {
        let specimen = mechanical.specimen();
        if !positive(temperature_k)
            || specimen.material().query().axes().get("T").copied() != Some(temperature_k)
        {
            return Err(admit("thermoelastic energy needs the candidate's actual source temperature"));
        }
        let StringPrestress::FixedThermalExtension { free_thermal_strain, .. } = specimen.prestress()
            else { return Err(admit("thermoelastic energy lost its fixed-support prescription")); };
        let expansion = specimen.thermal_expansion()
            .ok_or_else(|| admit("thermoelastic energy lost its original expansion path"))?;
        if expansion.selected_claim() != self.expansion.selected_claim()
            || expansion.reference().resolved().card_identity()
                != self.expansion.reference().resolved().card_identity()
            || expansion.reference().resolved().query_point()
                != self.expansion.reference().resolved().query_point()
            || expansion.reference().resolved().axis_quantities()
                != self.expansion.reference().resolved().axis_quantities()
            || expansion.current().resolved().card_identity() != specimen.material().card_identity()
            || expansion.current().resolved().query_point() != specimen.material().query_point()
            || expansion.current().resolved().axis_quantities() != specimen.material().axis_quantities()
            || expansion.free_linear_strain() != free_thermal_strain
        {
            return Err(admit("thermoelastic energy requires its unchanged original expansion source"));
        }
        let alpha = expansion.current().resolved()
            .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY)
            .ok_or_else(|| admit("thermoelastic energy lacks current source expansion"))?.value_si();
        let string = specimen.string();
        let modulus = specimen.material().property("young_modulus")
            .ok_or_else(|| admit("thermoelastic energy lacks Young's modulus"))?.value_si();
        let elastic = finite((string.length_m-self.reference_length_m)/self.reference_length_m
            -free_thermal_strain)?;
        let volume = finite(specimen.area_m2()*self.reference_length_m)?;
        let axial_free_energy_j = finite(0.5*volume*modulus*elastic*elastic)?;
        let mut entropy_j_k = finite(volume*(modulus*alpha*elastic
            -0.5*self.modulus_slope_pa_k*elastic*elastic))?;
        let mut modal_free_energy_j = 0.0;
        for (index, (mode, state)) in mechanical.modes().iter().zip(mechanical.states()).enumerate() {
            let k = finite((index+1) as f64*core::f64::consts::PI/string.length_m)?;
            let q = state.displacement_m_sqrt_kg;
            let derivative = finite((self.modulus_slope_pa_k*specimen.second_moment_m4()*k.powi(4)
                +specimen.area_m2()*(self.modulus_slope_pa_k*elastic-modulus*alpha)*k*k)
                /string.lin_density_kg_m)?;
            entropy_j_k = finite(entropy_j_k-0.5*derivative*q*q)?;
            let omega_q = finite(mode.angular_frequency_rad_s*q)?;
            modal_free_energy_j = finite(modal_free_energy_j+0.5*omega_q*omega_q)?;
        }
        Ok(ThermoelasticEnergy { axial_free_energy_j, modal_free_energy_j, entropy_j_k })
    }

    /// Compile bounds from the selected original table, once. The native
    /// source-line owner proves continuous support in the unchanged context.
    /// Eigenstrain and current alpha still come from their original owners.
    fn prepare_expansion_capacity(
        &self,
        cx: &fs_exec::Cx<'_>,
        mechanical: &LinearMaterialStringRuntime,
        card: &MaterialCard,
        curve: &EquilibriumEnthalpyPhaseCurve,
    ) -> Result<Option<Arc<[ExpansionCapacityBound]>>, ThermalStringError> {
        let claim = card.claims().claim(self.expansion.selected_claim())
            .ok_or_else(|| admit("retained expansion claim is absent"))?;
        let source = match &claim.value {
            PropertyValue::Scalar { .. } => return Ok(None),
            PropertyValue::Curve { abscissa, abscissa_dims, knots, .. }
                if abscissa == "T"
                    && *abscissa_dims == fs_qty::Temperature::DIMS
                    && claim.interpolation == fs_matdb::InterpolationPolicy::LinearInside
                    && (2..=MAX_HEAT_CAPACITY_ENTHALPY_KNOTS).contains(&knots.len()) => knots,
            _ => return Err(admit("live thermoelastic expansion must be scalar or a bounded continuous piecewise-affine temperature curve")),
        };
        let caloric = curve.knots();
        if !(2..=MAX_HEAT_CAPACITY_ENTHALPY_KNOTS).contains(&caloric.len()) {
            return Err(admit("thermoelastic expansion needs a bounded solid caloric interval"));
        }
        let lower = self.expansion.current().resolved().query()
            .with("T", caloric[0].temperature_k).map_err(ThermalStringError::Query)?;
        let upper = lower.clone().with("T", caloric[caloric.len()-1].temperature_k)
            .map_err(ThermalStringError::Query)?;
        card.claims().query_line_typed(&claim.key, "T", &lower, &upper,
            fs_matdb::ClaimSelection::Pinned(self.expansion.selected_claim()))
            .map_err(|_| admit("expansion source must continuously cover the original caloric interval"))?;
        if source.iter().all(|(_, alpha)| *alpha == self.expansion_per_k) {
            // An explicitly tabulated constant keeps the scalar capacity
            // admission, after its original bounded source support is checked.
            return Ok(None);
        }
        let specimen = mechanical.specimen();
        let volume = finite(specimen.area_m2()*self.reference_length_m)?;
        let modal_factor = finite(0.5*specimen.area_m2()/specimen.string().lin_density_kg_m)?;
        let b = self.modulus_slope_pa_k;
        let modulus = |temperature: f64| finite(
            b.mul_add(temperature-self.initial_temperature_k, self.initial_modulus_pa));
        let mut bounds = Vec::with_capacity(caloric.len()+source.len()-2);
        let mut source_index = 0;
        for (index, pair) in caloric.windows(2).enumerate() {
            check(cx)?;
            let mut lo = pair[0].temperature_k;
            let hi = pair[1].temperature_k;
            if !positive(hi-lo) {
                return Err(admit("thermoelastic caloric temperatures must strictly increase"));
            }
            let cp_minimum = curve.sensible_heat_capacity_bounds_j_kg_k(index)
                .map_err(|_| admit("string caloric heat capacity is unrepresentable"))?[0];
            while lo < hi {
                check(cx)?;
                while source_index+1 < source.len()-1 && source[source_index+1].0 <= lo {
                    source_index += 1;
                }
                let [(t0, a0), (t1, a1)] = [source[source_index], source[source_index+1]];
                let end = hi.min(t1);
                if lo < t0 || !positive(end-lo) || !positive(t1-t0) {
                    return Err(admit("expansion source does not cover a complete caloric segment"));
                }
                // The slope uses the ORIGINAL native segment, not differences
                // of clipped endpoint answers that can erase a small slope.
                let alpha_max = finite(a0.abs().max(a1.abs()))?;
                let slope_max = finite((a1-a0)/(t1-t0))?.abs();
                let e_max = finite(modulus(lo)?.abs().max(modulus(end)?.abs()))?;
                let cross = finite(2.0*b.abs()*alpha_max+e_max*slope_max)?;
                // Psi_TT = V(E alpha^2 - 2 E' alpha e - E alpha' e)
                //          - A/(2 mu) (2 E' alpha + E alpha') sum(k q)^2.
                // The native strain admission guarantees |e| <= strain_limit.
                // Endpoint extrema bound affine E and alpha; absolute values
                // retain both expansion and contraction without assuming signs.
                let outward = 1.0+256.0*f64::EPSILON;
                bounds.push(ExpansionCapacityBound {
                    upper_temperature_k: end,
                    caloric_capacity_j_k: finite(specimen.mass_kg()*cp_minimum)?,
                    axial_curvature_j_k2: finite(outward*volume
                        *(e_max*alpha_max*alpha_max+cross*self.strain_limit))?,
                    modal_curvature_per_gradient: finite(outward*modal_factor*cross)?,
                });
                lo = end;
            }
        }
        Ok(Some(bounds.into()))
    }

    /// Positive total fixed-q heat capacity over every original caloric segment.
    /// This checks the actual retained amplitudes, including E'*alpha and
    /// E*alpha' contributions. Scalar alpha keeps its exact quadratic minimum;
    /// varying alpha uses the prepared conservative strain-envelope bounds.
    /// It never claims tautness or mechanical source support outside an actual
    /// source-bound trial; those remain the native admissions.
    pub(crate) fn admit_capacity(
        &self,
        cx: &fs_exec::Cx<'_>,
        mechanical: &LinearMaterialStringRuntime,
        curve: &EquilibriumEnthalpyPhaseCurve,
    ) -> Result<(), ThermalStringError> {
        check(cx)?;
        let specimen = mechanical.specimen();
        let string = specimen.string();
        let knots = curve.knots();
        if knots.len() < 2 || knots.len() > MAX_HEAT_CAPACITY_ENTHALPY_KNOTS
            || knots.iter().any(|k| k.liquid_mass_fraction != 0.0
                || k.bulk_density_kg_m3 != knots[0].bulk_density_kg_m3)
        {
            return Err(admit("thermoelastic feedback requires a bounded single-solid constant-density caloric chart"));
        }
        let mut modal_gradient = 0.0;
        for (index, state) in mechanical.states().iter().enumerate() {
            let k = finite((index+1) as f64*core::f64::consts::PI/string.length_m)?;
            let kq = finite(k*state.displacement_m_sqrt_kg)?;
            modal_gradient = finite(modal_gradient+kq*kq)?;
        }
        if let Some(bounds) = &self.expansion_capacity {
            for bound in bounds.iter() {
                check(cx)?;
                let curvature = finite(bound.axial_curvature_j_k2
                    + bound.modal_curvature_per_gradient*modal_gradient)?;
                let capacity = finite(bound.caloric_capacity_j_k
                    - bound.upper_temperature_k*curvature)?;
                if !positive(capacity) {
                    return Err(admit("caloric capacity does not exceed the complete source-slope axial and modal thermoelastic bound"));
                }
            }
            return Ok(());
        }
        let volume = finite(specimen.area_m2()*self.reference_length_m)?;
        let alpha = self.expansion_per_k;
        let b = self.modulus_slope_pa_k;
        // Source Cp and Psi_TT are affine in T. Total capacity is quadratic;
        // endpoints and its stationary point retain its true minimum.
        let slope = finite(3.0*volume*b*alpha*alpha)?;
        let curvature = |temperature: f64| -> Result<f64, ThermalStringError> {
            let delta = finite(temperature-self.initial_temperature_k)?;
            let modulus = finite(b.mul_add(delta, self.initial_modulus_pa))?;
            let elastic = finite((-alpha).mul_add(delta, self.initial_elastic_strain))?;
            finite(volume*(modulus*alpha*alpha-2.0*b*alpha*elastic)
                -specimen.area_m2()/string.lin_density_kg_m*b*alpha*modal_gradient)
        };
        for (index, pair) in knots.windows(2).enumerate() {
            check(cx)?;
            let lo = pair[0].temperature_k;
            let hi = pair[1].temperature_k;
            let width = finite(hi-lo)?;
            if !positive(width) { return Err(admit("thermoelastic caloric temperatures must strictly increase")); }
            let bounds = curve.sensible_heat_capacity_bounds_j_kg_k(index)
                .map_err(|_|admit("string caloric heat capacity is unrepresentable"))?;
            let cp = curve.sensible_heat_capacities_on_interval(index)
                .unwrap_or([bounds[0]; 2]);
            let capacity_slope = finite(specimen.mass_kg()*((cp[1]-cp[0])/width))?;
            let total_capacity = |t: f64| -> Result<f64, ThermalStringError> {
                let fraction = (t-lo)/width;
                let cp_at_t = if fraction <= 0.5 {
                    (cp[1]-cp[0]).mul_add(fraction, cp[0])
                } else {
                    (cp[0]-cp[1]).mul_add((hi-t)/width, cp[1])
                };
                finite(specimen.mass_kg()*cp_at_t - t*curvature(t)?)
            };
            let lower_curvature = curvature(lo)?;
            let mut minimum = total_capacity(lo)?.min(total_capacity(hi)?);
            if slope != 0.0 {
                // An unrepresentably distant stationary point lies outside this
                // finite interval; the endpoint minimum remains sufficient.
                let stationary = 0.5*(lo+(capacity_slope-lower_curvature)/slope);
                if stationary.is_finite() && stationary > lo && stationary < hi {
                    minimum = minimum.min(total_capacity(stationary)?);
                }
            }
            if !positive(minimum) {
                return Err(admit("caloric capacity does not exceed the complete axial and modal thermoelastic contribution"));
            }
        }
        Ok(())
    }
}

fn admit(message: &'static str) -> ThermalStringError { ThermalStringError::Admission(message) }
fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn finite(value: f64) -> Result<f64, ThermalStringError> {
    if value.is_finite() { Ok(value) }
    else { Err(admit("unrepresentable thermoelastic string arithmetic")) }
}
fn check(cx: &fs_exec::Cx<'_>) -> Result<(), ThermalStringError> {
    cx.checkpoint().map_err(|_| ThermalStringError::Acoustic(super::AcousticRealizeError::Cancelled))
}
