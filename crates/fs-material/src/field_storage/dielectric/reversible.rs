//! Reversible isotropic dielectric response from one original monotone D(|E|).
//!
//! The voltage-form weak equation differentiates coenergy integral D dE. The
//! physical charge-form energy is integral E dD; neither equals E.D/2 for a
//! nonlinear law. Both are integrated directly on the same pinned source path.
//! Field and displacement inputs to evaluation are magnitudes. The signed
//! scalar work methods use the corresponding odd E(D) and even energy.
//!
//! These are fixed-context constitutive snapshots without remanence, memory,
//! leakage or a thermal/mechanical update. A source-dependent temperature
//! snapshot does not supply field entropy or internal-energy corrections.

use crate::{
    field_storage::FieldStorageError,
    state_point::{
        MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
        ScalarAdmissibility, ScalarPropertyRequirement, is_fixed_external_coordinate,
        resolve_material_state_point,
    },
};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyLineSlice, QueryPoint};
use fs_qty::{Dims, QuantitySpec};
use std::sync::Arc;

/// Electric-field magnitude in volts per metre.
pub const ELECTRIC_FIELD_DIMS: Dims = Dims([1, 1, -3, 0, -1, 0]);
/// Absolute displacement in coulombs per square metre, including vacuum response.
pub const ELECTRIC_DISPLACEMENT_DIMS: Dims = Dims([-2, 0, 1, 0, 1, 0]);
/// Explicit reversible, isotropic, zero-remanence displacement D(|E|).
/// Polarization alone, secant permittivity and hysteresis loops are not this key.
pub const REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY: &str = "reversible_electric_displacement";

/// One original source, context and finite field range, shared across cells.
/// Cloning shares immutable source receipts, slopes and exact segment integrals.
#[derive(Clone, Debug, PartialEq)]
pub struct ReversibleDielectricMaterial {
    prepared: Arc<Prepared>,
}

#[derive(Debug, PartialEq)]
struct Prepared {
    card: MaterialCard,
    material: ResolvedMaterialStatePoint,
    path: PropertyLineSlice,
    slopes_f_m: Vec<f64>,
    energy_at_knots_j_m3: Vec<f64>,
    coenergy_at_knots_j_m3: Vec<f64>,
}

/// Magnitude response and Legendre-dual potentials of the SAME source law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReversibleDielectricResponse {
    /// Nonnegative electric-field magnitude [V/m].
    pub field_strength_v_m: f64,
    /// Nonnegative absolute electric-displacement magnitude [C/m2].
    pub displacement_c_m2: f64,
    /// D/E [F/m], with the first segment's limit at zero.
    pub secant_permittivity_f_m: f64,
    /// dD/dE [F/m]; right derivative at interior knots, left at the last knot.
    pub differential_permittivity_f_m: f64,
    /// Physical fixed-charge storage integral_0^D E(s) ds [J/m3].
    pub energy_density_j_m3: f64,
    /// Voltage-form potential integral_0^E D(s) ds [J/m3].
    pub coenergy_density_j_m3: f64,
}

fn invalid(what: &'static str) -> FieldStorageError {
    FieldStorageError::Invalid(what)
}
fn finite(value: f64) -> Result<f64, FieldStorageError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "reversible dielectric arithmetic is unrepresentable",
        ))
    }
}
fn add(sum: &mut f64, correction: &mut f64, term: f64) -> Result<(), FieldStorageError> {
    let next = finite(*sum + term)?;
    *correction = finite(
        *correction
            + if sum.abs() >= term.abs() {
                (*sum - next) + term
            } else {
                (term - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}

impl ReversibleDielectricMaterial {
    /// Compile one continuously supported, source-selected D(E) curve or
    /// fixed-context surface slice from zero through the declared maximum.
    /// The existing line-query owner bounds the retained path to 4096 knots.
    ///
    /// The named field coordinate must already have explicit V/m dimensions.
    /// Every other coordinate must be a typed fixed external thermodynamic
    /// condition. D(0)=0 and positive representable dD/dE and dE/dD are required.
    /// A new context creates a new snapshot; it never changes live field history
    /// or accounts for parameter work or source-dependent thermal entropy.
    pub fn resolve(
        card: &MaterialCard,
        context: &QueryPoint,
        electric_field_axis: &str,
        maximum_field_v_m: f64,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FieldStorageError> {
        if context.axis_quantities().get(electric_field_axis)
            != Some(&QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS))
            || !context.axes().contains_key(electric_field_axis)
            || !maximum_field_v_m.is_finite()
            || maximum_field_v_m <= 0.0
        {
            return Err(invalid(
                "reversible dielectric needs a typed V/m coordinate and positive finite range",
            ));
        }
        let point = |field| {
            context
                .clone()
                .with(electric_field_axis, field)
                .map_err(|_| invalid("dielectric field query is unrepresentable"))
        };
        let lower = point(0.0)?;
        let upper = point(maximum_field_v_m)?;
        let external = |axis: &str| {
            lower.axis_quantities().contains_key(axis)
                && is_fixed_external_coordinate(
                    axis,
                    lower.axis_quantities().get(axis).copied(),
                    lower.axes().get(axis).copied(),
                )
        };
        if lower
            .axes()
            .keys()
            .any(|axis| axis != electric_field_axis && !external(axis))
        {
            return Err(invalid(
                "dielectric context contains an unowned response coordinate",
            ));
        }
        let requirement = ScalarPropertyRequirement::try_new(
            REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
            ELECTRIC_DISPLACEMENT_DIMS,
            ScalarAdmissibility::NonNegative,
        )?;
        let material = resolve_material_state_point(card, &lower, &[requirement], selection)?;
        let property = material
            .property(REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY)
            .ok_or(invalid("reversible displacement source is missing"))?;
        let selected = property.answer().receipt.selected;
        let claim = card
            .claims()
            .claim(selected)
            .ok_or(invalid("original reversible displacement claim is absent"))?;
        if !property
            .source_axes()
            .iter()
            .any(|(axis, dims)| axis == electric_field_axis && *dims == ELECTRIC_FIELD_DIMS)
            || property.source_axes().iter().any(|(axis, dims)| {
                if axis == electric_field_axis {
                    *dims != ELECTRIC_FIELD_DIMS
                } else {
                    !external(axis)
                        || lower
                            .axis_quantities()
                            .get(axis)
                            .is_none_or(|quantity| quantity.dims() != *dims)
                }
            })
            || claim
                .validity
                .bounds()
                .keys()
                .any(|axis| axis != electric_field_axis && !external(axis))
        {
            return Err(invalid(
                "reversible dielectric source requires one owned E axis and fixed external conditions",
            ));
        }
        let path = card
            .claims()
            .query_line_typed(
                &claim.key,
                electric_field_axis,
                &lower,
                &upper,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY.to_owned(),
                source: Box::new(source),
            })?;
        if path.knots.len() < 2
            || path.knots[0].0 != 0.0
            || path.knots[0].1.evidence.value.value != 0.0
            || path.knots.last().map(|k| k.0) != Some(maximum_field_v_m)
        {
            return Err(invalid(
                "reversible dielectric must include D(0)=0 and the complete source range",
            ));
        }
        let mut slopes_f_m = Vec::with_capacity(path.knots.len() - 1);
        let mut energy_at_knots_j_m3 = vec![0.0];
        let mut coenergy_at_knots_j_m3 = vec![0.0];
        let (mut energy, mut energy_correction) = (0.0, 0.0);
        let (mut coenergy, mut coenergy_correction) = (0.0, 0.0);
        for pair in path.knots.windows(2) {
            let (e0, e1) = (pair[0].0, pair[1].0);
            let (d0, d1) = (
                pair[0].1.evidence.value.value,
                pair[1].1.evidence.value.value,
            );
            let (de, dd) = (e1 - e0, d1 - d0);
            let slope = dd / de;
            let energy_increment = f64::midpoint(e0, e1) * dd;
            let coenergy_increment = f64::midpoint(d0, d1) * de;
            if [
                de,
                dd,
                slope,
                slope.recip(),
                energy_increment,
                coenergy_increment,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            {
                return Err(invalid(
                    "dielectric source requires positive representable differential permittivity and both energies",
                ));
            }
            add(&mut energy, &mut energy_correction, energy_increment)?;
            add(&mut coenergy, &mut coenergy_correction, coenergy_increment)?;
            slopes_f_m.push(slope);
            energy_at_knots_j_m3.push(finite(energy + energy_correction)?);
            coenergy_at_knots_j_m3.push(finite(coenergy + coenergy_correction)?);
        }
        Ok(Self {
            prepared: Arc::new(Prepared {
                card: card.clone(),
                material,
                path,
                slopes_f_m,
                energy_at_knots_j_m3,
                coenergy_at_knots_j_m3,
            }),
        })
    }

    /// Original immutable card; no fitted or synthesized source is substituted.
    pub fn card(&self) -> &MaterialCard {
        &self.prepared.card
    }
    /// Original zero-field snapshot with the selected claim and complete context.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.prepared.material
    }
    /// Every canonical source breakpoint and its original answer/receipt.
    pub fn source_path(&self) -> &PropertyLineSlice {
        &self.prepared.path
    }
    /// Largest admitted field magnitude [V/m], never extrapolated.
    pub fn maximum_field_strength_v_m(&self) -> f64 {
        self.prepared.path.knots.last().unwrap().0
    }
    /// Largest admitted displacement magnitude [C/m2] on the same source path.
    pub fn maximum_displacement_c_m2(&self) -> f64 {
        self.prepared
            .path
            .knots
            .last()
            .unwrap()
            .1
            .evidence
            .value
            .value
    }
    fn admit_displacement(&self, displacement: f64) -> Result<(), FieldStorageError> {
        if displacement.is_finite() && displacement.abs() <= self.maximum_displacement_c_m2() {
            Ok(())
        } else {
            Err(invalid(
                "electric displacement leaves the original reversible source range",
            ))
        }
    }
    fn displacement_segment(&self, displacement: f64) -> usize {
        self.prepared
            .path
            .knots
            .partition_point(|(_, answer)| answer.evidence.value.value <= displacement)
            .saturating_sub(1)
            .min(self.prepared.slopes_f_m.len() - 1)
    }

    /// Evaluate a NONNEGATIVE field magnitude with no query or allocation.
    /// O(log knots) binary search selects the original piecewise-affine segment.
    pub fn evaluate(
        &self,
        field_v_m: f64,
    ) -> Result<ReversibleDielectricResponse, FieldStorageError> {
        if !field_v_m.is_finite()
            || field_v_m < 0.0
            || field_v_m > self.maximum_field_strength_v_m()
        {
            return Err(invalid(
                "electric field leaves the original reversible source range",
            ));
        }
        let i = self
            .prepared
            .path
            .knots
            .partition_point(|(e, _)| *e <= field_v_m)
            .saturating_sub(1)
            .min(self.prepared.slopes_f_m.len() - 1);
        let (e0, answer) = &self.prepared.path.knots[i];
        let d0 = answer.evidence.value.value;
        let de = field_v_m - e0;
        let dd = self.prepared.slopes_f_m[i] * de;
        let displacement = if field_v_m == self.maximum_field_strength_v_m() {
            self.maximum_displacement_c_m2()
        } else {
            self.prepared.slopes_f_m[i].mul_add(de, d0)
        };
        self.response(i, field_v_m, displacement, de, dd)
    }

    /// Invert and evaluate at a NONNEGATIVE displacement magnitude. The SAME
    /// source segment supplies E(D), dE/dD and both potentials; no inverse fit.
    pub fn evaluate_displacement(
        &self,
        displacement_c_m2: f64,
    ) -> Result<ReversibleDielectricResponse, FieldStorageError> {
        self.admit_displacement(displacement_c_m2)?;
        if displacement_c_m2 < 0.0 {
            return Err(invalid(
                "dielectric evaluation requires a displacement magnitude",
            ));
        }
        let i = self.displacement_segment(displacement_c_m2);
        let (e0, answer) = &self.prepared.path.knots[i];
        let dd = displacement_c_m2 - answer.evidence.value.value;
        let de = dd / self.prepared.slopes_f_m[i];
        let e1 = self.prepared.path.knots[i + 1].0;
        let field = if displacement_c_m2 == self.maximum_displacement_c_m2() {
            self.maximum_field_strength_v_m()
        } else {
            // Only roundoff of an already admitted affine inverse is limited.
            (e0 + de).clamp(*e0, e1)
        };
        self.response(i, field, displacement_c_m2, de, dd)
    }

    /// Nonnegative E(|D|) [V/m], inverting the original admitted source exactly.
    pub fn field_at_displacement(&self, displacement_c_m2: f64) -> Result<f64, FieldStorageError> {
        Ok(self
            .evaluate_displacement(displacement_c_m2)?
            .field_strength_v_m)
    }

    fn response(
        &self,
        i: usize,
        e: f64,
        d: f64,
        de: f64,
        dd: f64,
    ) -> Result<ReversibleDielectricResponse, FieldStorageError> {
        let (e0, answer) = &self.prepared.path.knots[i];
        let d0 = answer.evidence.value.value;
        // Work uses de and dd, not subtraction of large rounded total fields.
        let (energy, coenergy) =
            if e == self.maximum_field_strength_v_m() && d == self.maximum_displacement_c_m2() {
                (
                    *self.prepared.energy_at_knots_j_m3.last().unwrap(),
                    *self.prepared.coenergy_at_knots_j_m3.last().unwrap(),
                )
            } else {
                (
                    (e0 + 0.5 * de).mul_add(dd, self.prepared.energy_at_knots_j_m3[i]),
                    (d0 + 0.5 * dd).mul_add(de, self.prepared.coenergy_at_knots_j_m3[i]),
                )
            };
        let differential = self.prepared.slopes_f_m[i];
        let secant = if e == 0.0 { differential } else { d / e };
        if [e, d, secant, energy, coenergy]
            .iter()
            .any(|v| !v.is_finite())
            || secant <= 0.0
            || (e > 0.0 || d > 0.0) && (e <= 0.0 || d <= 0.0 || energy <= 0.0 || coenergy <= 0.0)
        {
            return Err(invalid("dielectric response or energy is unrepresentable"));
        }
        Ok(ReversibleDielectricResponse {
            field_strength_v_m: e,
            displacement_c_m2: d,
            secant_permittivity_f_m: secant,
            differential_permittivity_f_m: differential,
            energy_density_j_m3: energy,
            coenergy_density_j_m3: coenergy,
        })
    }

    // Walk the odd inverse law in the requested direction, splitting at zero
    // and original signed displacement knots. Each interval has dE/dD > 0.
    fn walk_displacement(
        &self,
        before: f64,
        after: f64,
        mut visit: impl FnMut(f64, f64) -> Result<(), FieldStorageError>,
    ) -> Result<(), FieldStorageError> {
        let direction = if after >= before { 1.0 } else { -1.0 };
        let mut position = before;
        while position != after {
            let magnitude = position.abs();
            let outward = position == 0.0 || position.signum() == direction;
            let mut i = self.displacement_segment(magnitude);
            if !outward && i > 0 && magnitude == self.prepared.path.knots[i].1.evidence.value.value
            {
                i -= 1;
            }
            let boundary_magnitude = self.prepared.path.knots[i + usize::from(outward)]
                .1
                .evidence
                .value
                .value;
            let boundary = if position == 0.0 {
                direction
            } else {
                position.signum()
            } * boundary_magnitude;
            let next = if direction > 0.0 {
                after.min(boundary)
            } else {
                after.max(boundary)
            };
            visit(next - position, self.prepared.slopes_f_m[i].recip())?;
            position = next;
        }
        Ok(())
    }

    fn integrate_relative_field(
        &self,
        before: f64,
        after: f64,
        initial_field: f64,
    ) -> Result<f64, FieldStorageError> {
        let (mut field, mut field_correction) = (initial_field, 0.0);
        let (mut work, mut work_correction) = (0.0, 0.0);
        self.walk_displacement(before, after, |dd, inverse| {
            let de = finite(dd * inverse)?;
            let midpoint = finite(field + finite(field_correction + 0.5 * de)?)?;
            let term = finite(dd * midpoint)?;
            if dd != 0.0 && (de == 0.0 || midpoint != 0.0 && term == 0.0) {
                return Err(invalid("nonzero dielectric work increment underflows"));
            }
            add(&mut work, &mut work_correction, term)?;
            add(&mut field, &mut field_correction, de)
        })?;
        finite(work + work_correction)
    }

    /// Stable w(after)-w(before) [J/m3] for SIGNED scalar displacement. Energy
    /// is even. Only the interval between magnitudes is integrated, so reversal
    /// at equal magnitude returns exactly zero without cancelling stored bias.
    pub fn energy_change_density_j_m3(
        &self,
        before_d: f64,
        after_d: f64,
    ) -> Result<f64, FieldStorageError> {
        self.admit_displacement(before_d)?;
        self.admit_displacement(after_d)?;
        let before = before_d.abs();
        let after = after_d.abs();
        if before == after {
            return Ok(0.0);
        }
        let field = self.field_at_displacement(before)?;
        self.integrate_relative_field(before, after, field)
    }

    /// Stable coenergy(after)-coenergy(before) [J/m3] for NONNEGATIVE electric
    /// field magnitudes. Original crossed-segment D dE integrals preserve small
    /// Newton line-search work beneath a large unchanged field-energy bias.
    pub fn coenergy_change_density_j_m3(
        &self,
        before_e: f64,
        after_e: f64,
    ) -> Result<f64, FieldStorageError> {
        for e in [before_e, after_e] {
            if !e.is_finite() || e < 0.0 || e > self.maximum_field_strength_v_m() {
                return Err(invalid(
                    "electric field leaves the original reversible source range",
                ));
            }
        }
        let increasing = after_e >= before_e;
        let mut position = before_e;
        let (mut work, mut correction) = (0.0, 0.0);
        while position != after_e {
            let mut i = self
                .prepared
                .path
                .knots
                .partition_point(|(e, _)| *e <= position)
                .saturating_sub(1)
                .min(self.prepared.slopes_f_m.len() - 1);
            if !increasing && i > 0 && position == self.prepared.path.knots[i].0 {
                i -= 1;
            }
            let (e0, answer) = &self.prepared.path.knots[i];
            let next = if increasing {
                after_e.min(self.prepared.path.knots[i + 1].0)
            } else {
                after_e.max(*e0)
            };
            let de = next - position;
            let midpoint_d = finite(
                self.prepared.slopes_f_m[i]
                    .mul_add((position - e0) + 0.5 * de, answer.evidence.value.value),
            )?;
            let term = finite(de * midpoint_d)?;
            if term == 0.0 {
                return Err(invalid("nonzero dielectric coenergy change underflows"));
            }
            add(&mut work, &mut correction, term)?;
            position = next;
        }
        finite(work + correction)
    }

    /// Stable isotropic coenergy change between actual field VECTORS [V/m].
    ///
    /// The source integral supplies the mean radial permittivity on the same
    /// crossed segments. It weights the componentwise change of |E|²/2, keeping
    /// transverse work even when rounded endpoint magnitudes coincide. The
    /// local D/E limit is used in that case. Corrected midpoint products apply
    /// also when the endpoint magnitudes differ, rather than relying on their
    /// rounded squared-norm difference. All reductions are allocation-free.
    ///
    /// Multiplication by permittivity precedes squared products, so a raw field
    /// square need not be representable when the physical work is. Magnitudes
    /// still obey the original source support and floating-point admission;
    /// this supplies no interval enclosure of a field near a source boundary.
    pub fn coenergy_change_vector_density_j_m3<const N: usize>(
        &self,
        before: [f64; N],
        after: [f64; N],
    ) -> Result<f64, FieldStorageError> {
        if before.iter().chain(&after).any(|v| !v.is_finite()) {
            return Err(invalid("dielectric field vectors must be finite"));
        }
        let before_magnitude = before.iter().fold(0.0_f64, |m, v| m.hypot(*v));
        let after_magnitude = after.iter().fold(0.0_f64, |m, v| m.hypot(*v));
        let radial_permittivity = if before_magnitude == after_magnitude {
            self.evaluate(before_magnitude)?.secant_permittivity_f_m
        } else {
            let work = self.coenergy_change_density_j_m3(before_magnitude, after_magnitude)?;
            // Do not form (after²-before²)/2: either square can overflow,
            // underflow, or erase the small increment that this method retains.
            let mean_displacement = finite(work / (after_magnitude - before_magnitude))?;
            finite(mean_displacement / f64::midpoint(before_magnitude, after_magnitude))?
        };
        if radial_permittivity <= 0.0 {
            return Err(invalid(
                "radial dielectric work coefficient is unrepresentable",
            ));
        }
        let (mut work, mut correction) = (0.0, 0.0);
        for (a, b) in before.into_iter().zip(after) {
            // Halves keep both sum and difference finite even for opposite
            // near-maximum components. TwoSum preserves each midpoint error.
            let (midpoint, midpoint_error) = two_sum(0.5 * a, 0.5 * b);
            let (half_delta, half_delta_error) = two_sum(0.5 * b, -0.5 * a);
            for delta in [half_delta, half_delta_error] {
                for middle in [midpoint, midpoint_error] {
                    add_twice_weighted_product(
                        &mut work,
                        &mut correction,
                        delta,
                        middle,
                        radial_permittivity,
                    )?;
                }
            }
        }
        finite(work + correction)
    }

    /// Stable integral_before^after [E(s)-E(reference)] ds [J/m3], accepting
    /// SIGNED scalar displacement on the odd inverse law. This supplies native
    /// discrete gradients and midpoint remainders without subtracting large
    /// endpoint energies or nearly equal biased electric fields. The result is
    /// signed in general and is not physical dielectric dissipation.
    pub fn energy_change_remainder_density_j_m3(
        &self,
        before_d: f64,
        after_d: f64,
        reference_d: f64,
    ) -> Result<f64, FieldStorageError> {
        for d in [before_d, after_d, reference_d] {
            self.admit_displacement(d)?;
        }
        let (mut offset, mut correction) = (0.0, 0.0);
        self.walk_displacement(reference_d, before_d, |dd, inverse| {
            let increment = finite(dd * inverse)?;
            if dd != 0.0 && increment == 0.0 {
                return Err(invalid("nonzero dielectric field difference underflows"));
            }
            add(&mut offset, &mut correction, increment)
        })?;
        self.integrate_relative_field(before_d, after_d, finite(offset + correction)?)
    }
}

// Inputs are half-components; their sum cannot overflow for finite vectors.
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let sum = a + b;
    let recovered_b = sum - a;
    (sum, (a - (sum - recovered_b)) + (b - recovered_b))
}

fn add_twice_weighted_product(
    sum: &mut f64,
    correction: &mut f64,
    mut a: f64,
    mut b: f64,
    permittivity: f64,
) -> Result<(), FieldStorageError> {
    if a == 0.0 || b == 0.0 {
        return Ok(());
    }
    // For a half-difference and midpoint at least one magnitude can be doubled
    // without overflow. Double the smaller, retaining subnormal work too.
    if a.abs() <= b.abs() {
        a = finite(2.0 * a)?;
    } else {
        b = finite(2.0 * b)?;
    }
    // Small permittivity first scales the larger field; large permittivity
    // first scales the smaller. This avoids unnecessary intermediate range loss.
    let (scaled, other) = if (a.abs() >= b.abs()) == (permittivity < 1.0) {
        (a, b)
    } else {
        (b, a)
    };
    let displacement = finite(permittivity * scaled)?;
    let displacement_error = finite(permittivity.mul_add(scaled, -displacement))?;
    let product = finite(displacement * other)?;
    if product == 0.0 {
        return Err(invalid(
            "nonzero vector dielectric work increment underflows",
        ));
    }
    add(sum, correction, product)?;
    add(
        sum,
        correction,
        finite(displacement.mul_add(other, -product))?,
    )?;
    add(sum, correction, finite(displacement_error * other)?)
}
