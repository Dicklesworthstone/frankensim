//! Evaluated fission energy on the original reaction-weighted condensation path.
//!
//! A source reference point does not specify its energy dependence. Callers
//! explicitly reject missing laws, freeze the reference, or select the ENDF-102
//! (2023) Sher--Beck prescription with an independently supplied nu(0). Tables
//! override that prescription component by component. Source ER/ET remain
//! independent observations; they are not silently replaced by component sums.
//!
//! Tables, linear polynomials and Sher--Beck segments are retained exactly up to
//! arithmetic. Higher polynomials are linearized under the interpolation bound
//! h^2 max|q''| / 8. The reported bound concerns that approximation in exact
//! arithmetic, not nuclear-data uncertainty or a certified floating-point bound.
//! Positive reaction weighting carries this uniform error bound into each
//! group's mean. No midpoint energy replaces the reaction integral.

use std::sync::Arc;

use super::{FissionEnergyComponent, FissionEnergyRelease};
use crate::cross_sections::NeutronReaction;
use crate::endf::neutron_yield::{FissionNeutronKind, FissionNeutronYield};
use crate::endf::reaction_weighted::{
    ReactionWeightedCollapse, ReactionWeightedLimits, TabulatedYield, TabulatedYieldData,
    TabulatedYieldLimits, collapse_reaction_weighted,
};
use crate::endf::spectrum::FluxSpectrum;
use crate::endf::{
    EV_J, EndfError, EndfInterpolation, EndfInterpolationRegion, EndfMetadata,
    PointwiseCrossSection, budget, convex,
};
use crate::irradiation::NeutronGroups;
use crate::{FissionError, nonnegative};

/// Explicit continuation for NPLY=0 components lacking a TAB1 law.
#[derive(Clone, Debug)]
pub enum ReferenceEnergyTreatment {
    /// Every component must supply a polynomial or a tabulated law.
    Reject,
    /// Caller-declared constant approximation over the requested domain.
    ConstantReference,
    /// ENDF-102 (2023) section 1.5.2.1. The supplied prompt or total nu curve
    /// must belong to the same evaluation and cover the entire domain. nu(0)
    /// is an explicit input; an unavailable zero-energy value is not extrapolated.
    SherBeck2023 {
        neutron_yield: Arc<FissionNeutronYield>,
        reference_nu_at_zero: f64,
    },
}

/// Shared preparation bounds across all nine components and the incident curve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FissionEnergyPreparationLimits {
    pub max_curve_points: usize,
    pub max_evaluations: usize,
    pub max_depth: usize,
    /// Allowed polynomial interpolation error in J per fission.
    pub absolute_tolerance_j: f64,
    /// Multiplies the smaller interval endpoint energy, not a signed coefficient.
    pub relative_tolerance: f64,
}

impl Default for FissionEnergyPreparationLimits {
    fn default() -> Self {
        Self {
            max_curve_points: 1_000_000,
            max_evaluations: 4_000_000,
            max_depth: 32,
            absolute_tolerance_j: 1e-6 * EV_J,
            relative_tolerance: 1e-8,
        }
    }
}

/// Shared immutable component curves, normalized by the explicit energy 1 eV.
#[derive(Clone, Debug)]
pub struct PreparedFissionEnergy {
    evaluation: Arc<FissionEnergyRelease>,
    domain_j: [f64; 2],
    reference: ReferenceEnergyTreatment,
    curves: Arc<[TabulatedYield]>,
    interpolation_error_j: [f64; 9],
}

/// Group means of each original component, with the independently weighted
/// incident-neutron energy. None denotes a group with no fission reactions.
#[derive(Clone, Debug)]
pub struct CollapsedFissionEnergy {
    evaluation: Arc<FissionEnergyRelease>,
    source_preparation: PreparedFissionEnergy,
    weighted: ReactionWeightedCollapse,
    source_cross_section: PointwiseCrossSection,
    source_spectrum: FluxSpectrum,
    components_j: [Vec<Option<f64>>; 9],
    incident_energy_j: Vec<Option<f64>>,
    interpolation_error_j: [f64; 9],
}

impl PreparedFissionEnergy {
    pub fn new(
        evaluation: Arc<FissionEnergyRelease>,
        domain_j: [f64; 2],
        reference: ReferenceEnergyTreatment,
        limits: FissionEnergyPreparationLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        if !domain_j.into_iter().all(nonnegative)
            || domain_j[0] >= domain_j[1]
            || !nonnegative(limits.absolute_tolerance_j)
            || !nonnegative(limits.relative_tolerance)
        {
            return Err(invalid(
                "finite increasing energy domain and nonnegative interpolation allowances required",
            ));
        }
        if let ReferenceEnergyTreatment::SherBeck2023 {
            neutron_yield,
            reference_nu_at_zero,
        } = &reference
        {
            if !same_evaluation(evaluation.metadata(), neutron_yield.metadata())
                || neutron_yield.kind() == FissionNeutronKind::Delayed
                || !nonnegative(*reference_nu_at_zero)
            {
                return Err(invalid(
                    "Sher-Beck requires same-evaluation prompt/total nu and explicit finite nu(0)",
                ));
            }
            neutron_yield.curve().value_at_energy(domain_j[0])?;
            neutron_yield.curve().value_at_energy(domain_j[1])?;
        }
        let mut work = Work {
            points: 0,
            evaluations: 0,
            limits,
        };
        let mut curves = Vec::with_capacity(10);
        let mut interpolation_error_j = [0.0; 9];
        for component in FissionEnergyComponent::ALL {
            checkpoint()?;
            let curve = if let Some(table) = evaluation.tabulation(component) {
                table.value_at_energy_j(domain_j[0])?;
                table.value_at_energy_j(domain_j[1])?;
                work.points(table.curve.energies_j().len())?;
                table.curve.clone()
            } else if evaluation.polynomial_order() > 0 {
                let (curve, error) =
                    polynomial_curve(&evaluation, component, domain_j, &mut work, &mut checkpoint)?;
                interpolation_error_j[component.index()] = error;
                curve
            } else {
                match &reference {
                    ReferenceEnergyTreatment::Reject => {
                        return Err(invalid(
                            "a reference-only fission component needs an explicit continuation policy",
                        ));
                    }
                    ReferenceEnergyTreatment::ConstantReference => {
                        work.points(2)?;
                        make_curve(
                            &evaluation,
                            component,
                            domain_j.to_vec(),
                            vec![evaluation.reference_value_j(component) / EV_J; 2],
                            &mut checkpoint,
                        )?
                    }
                    ReferenceEnergyTreatment::SherBeck2023 {
                        neutron_yield,
                        reference_nu_at_zero,
                    } => sher_beck_curve(
                        &evaluation,
                        component,
                        domain_j,
                        neutron_yield,
                        *reference_nu_at_zero,
                        &mut work,
                        &mut checkpoint,
                    )?,
                }
            };
            curves.push(curve);
        }
        // This independent curve preserves the incident kinetic-energy mean.
        // Outgoing product energy and net nuclear Q must not be conflated.
        work.points(2)?;
        curves.push(TabulatedYield::new_with_checkpoint(
            TabulatedYieldData {
                data_id: format!("{} | incident energy / 1 eV", evaluation.data_id()),
                energies_j: domain_j.to_vec(),
                values: domain_j.iter().map(|value| value / EV_J).collect(),
                interpolation_regions: vec![EndfInterpolationRegion {
                    end_point: 2,
                    interpolation: EndfInterpolation::Linear,
                }],
            },
            TabulatedYieldLimits {
                max_points: 2,
                max_interpolation_regions: 1,
            },
            &mut checkpoint,
        )?);
        checkpoint()?;
        Ok(Self {
            evaluation,
            domain_j,
            reference,
            curves: curves.into(),
            interpolation_error_j,
        })
    }

    pub fn evaluation(&self) -> &Arc<FissionEnergyRelease> {
        &self.evaluation
    }
    pub fn domain_j(&self) -> [f64; 2] {
        self.domain_j
    }
    pub fn reference_treatment(&self) -> &ReferenceEnergyTreatment {
        &self.reference
    }
    /// The admitted energy/1-eV laws, shared with each collapsed temperature.
    pub(super) fn component_curves(&self) -> &[TabulatedYield] {
        &self.curves
    }
    pub fn interpolation_error_j(&self, component: FissionEnergyComponent) -> f64 {
        self.interpolation_error_j[component.index()]
    }

    /// Evaluate the admitted prepared approximation without extrapolation.
    pub fn value_at_energy_j(
        &self,
        component: FissionEnergyComponent,
        energy_j: f64,
    ) -> Result<f64, EndfError> {
        self.check_domain(energy_j, energy_j)?;
        to_joules(self.curves[component.index()].value_at_energy(energy_j)?)
    }

    /// Preserve integral(sigma_f * phi * component energy) with the existing
    /// positive knot-union integrator. Flux amplitudes may later vary between
    /// groups; the within-group spectrum and processed state remain fixed.
    pub fn collapse(
        &self,
        fission: &PointwiseCrossSection,
        spectrum: &FluxSpectrum,
        groups: Arc<NeutronGroups>,
        limits: ReactionWeightedLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<CollapsedFissionEnergy, EndfError> {
        checkpoint()?;
        if fission.reaction() != NeutronReaction::Fission
            || !same_evaluation(self.evaluation.metadata(), fission.metadata())
        {
            return Err(invalid(
                "fission energy and cross section require the same isotope, target mass and evaluation revision",
            ));
        }
        let edges = groups.edges_j();
        self.check_domain(edges[0], edges[edges.len() - 1])?;
        let weighted = collapse_reaction_weighted(
            fission,
            spectrum,
            &self.curves,
            groups,
            limits,
            &mut checkpoint,
        )?;
        let mut components_j: [Vec<Option<f64>>; 9] = std::array::from_fn(|_| Vec::new());
        for component in FissionEnergyComponent::ALL {
            let row = &weighted.yields[component.index()].values;
            let values = &mut components_j[component.index()];
            values.reserve(row.len());
            for &value in row {
                checkpoint()?;
                values.push(value.map(to_joules).transpose()?);
            }
        }
        let mut incident_energy_j = Vec::with_capacity(weighted.flux.groups().len());
        for &value in &weighted.yields[9].values {
            checkpoint()?;
            incident_energy_j.push(value.map(to_joules).transpose()?);
        }
        checkpoint()?;
        Ok(CollapsedFissionEnergy {
            evaluation: Arc::clone(&self.evaluation),
            source_preparation: self.clone(),
            weighted,
            components_j,
            source_cross_section: fission.clone(),
            source_spectrum: spectrum.clone(),
            incident_energy_j,
            interpolation_error_j: self.interpolation_error_j,
        })
    }

    fn check_domain(&self, lower: f64, upper: f64) -> Result<(), EndfError> {
        if !nonnegative(lower)
            || !nonnegative(upper)
            || lower < self.domain_j[0]
            || upper > self.domain_j[1]
        {
            Err(invalid("energy outside the prepared fission-energy domain"))
        } else {
            Ok(())
        }
    }
}

impl CollapsedFissionEnergy {
    /// Exact prepared component laws and reference treatment used for collapse.
    /// This prevents a temperature series from silently changing energy policy.
    pub fn source_preparation(&self) -> &PreparedFissionEnergy {
        &self.source_preparation
    }

    pub fn evaluation(&self) -> &Arc<FissionEnergyRelease> {
        &self.evaluation
    }
    /// Original pointwise reaction, retained by shared immutable arrays.
    pub fn source_cross_section(&self) -> &PointwiseCrossSection {
        &self.source_cross_section
    }
    /// Original within-group spectral shape, not merely its integrated flux.
    pub fn source_spectrum(&self) -> &FluxSpectrum {
        &self.source_spectrum
    }
    /// Numerical ordinates in this original result are energy divided by 1 eV.
    /// Use component_j for physical SI means, including independent ER and ET.
    pub fn reaction_weighted(&self) -> &ReactionWeightedCollapse {
        &self.weighted
    }
    pub fn component_j(&self, component: FissionEnergyComponent) -> &[Option<f64>] {
        &self.components_j[component.index()]
    }
    pub fn incident_energy_j(&self) -> &[Option<f64>] {
        &self.incident_energy_j
    }
    pub fn interpolation_error_j(&self, component: FissionEnergyComponent) -> f64 {
        self.interpolation_error_j[component.index()]
    }
}

fn same_evaluation(a: &EndfMetadata, b: &EndfMetadata) -> bool {
    a.nuclide == b.nuclide
        && a.atomic_weight_ratio == b.atomic_weight_ratio
        && a.library_id == b.library_id
        && a.library_version == b.library_version
        && a.library_release == b.library_release
        && a.material_revision == b.material_revision
}

struct Work {
    points: usize,
    evaluations: usize,
    limits: FissionEnergyPreparationLimits,
}
impl Work {
    fn points(&mut self, count: usize) -> Result<(), EndfError> {
        self.points = self
            .points
            .checked_add(count)
            .filter(|&n| n <= self.limits.max_curve_points)
            .ok_or_else(|| budget("fission energy curve points"))?;
        Ok(())
    }
    fn evaluate(
        &mut self,
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(), EndfError> {
        checkpoint()?;
        self.evaluations = self
            .evaluations
            .checked_add(1)
            .filter(|&n| n <= self.limits.max_evaluations)
            .ok_or_else(|| budget("fission energy evaluations"))?;
        Ok(())
    }
}

fn make_curve(
    evaluation: &FissionEnergyRelease,
    component: FissionEnergyComponent,
    energies_j: Vec<f64>,
    values_ev: Vec<f64>,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<TabulatedYield, EndfError> {
    let count = energies_j.len();
    TabulatedYield::new_with_checkpoint(
        TabulatedYieldData {
            data_id: format!(
                "{} | MF1/MT458 {:?} / 1 eV",
                evaluation.data_id(),
                component
            ),
            energies_j,
            values: values_ev,
            interpolation_regions: vec![EndfInterpolationRegion {
                end_point: count,
                interpolation: EndfInterpolation::Linear,
            }],
        },
        TabulatedYieldLimits {
            max_points: count,
            max_interpolation_regions: 1,
        },
        checkpoint,
    )
}

fn polynomial_curve(
    evaluation: &FissionEnergyRelease,
    component: FissionEnergyComponent,
    domain: [f64; 2],
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(TabulatedYield, f64), EndfError> {
    work.evaluate(checkpoint)?;
    let ya = evaluation.value_at_energy_j(component, domain[0])?;
    work.evaluate(checkpoint)?;
    let yb = evaluation.value_at_energy_j(component, domain[1])?;
    work.points(1)?;
    let mut energies = vec![domain[0]];
    let mut values = vec![ya / EV_J];
    let positive_coefficients = evaluation
        .coefficients(component)
        .iter()
        .all(|c| c.value >= 0.0);
    let mut intervals = vec![(domain[0], domain[1], ya, yb, 0_usize)];
    let mut maximum_error = 0.0_f64;
    while let Some((a, b, left, right, depth)) = intervals.pop() {
        checkpoint()?;
        let error = polynomial_error(evaluation, component, a, b)?;
        let tolerance =
            work.limits.absolute_tolerance_j + work.limits.relative_tolerance * left.min(right);
        if !tolerance.is_finite() {
            return Err(FissionError::Arithmetic.into());
        }
        if error <= tolerance && (positive_coefficients || left.min(right) >= error) {
            work.points(1)?;
            energies.push(b);
            values.push(right / EV_J);
            maximum_error = maximum_error.max(error);
            continue;
        }
        if depth >= work.limits.max_depth {
            return Err(budget(
                "fission-energy polynomial refinement depth or nonnegative domain",
            ));
        }
        let middle = a + (b - a) * 0.5;
        if middle <= a || middle >= b {
            return Err(FissionError::Arithmetic.into());
        }
        work.evaluate(checkpoint)?;
        let value = evaluation.value_at_energy_j(component, middle)?;
        intervals.push((middle, b, value, right, depth + 1));
        intervals.push((a, middle, left, value, depth + 1));
    }
    Ok((
        make_curve(evaluation, component, energies, values, checkpoint)?,
        maximum_error,
    ))
}

// For nonnegative incident energies, the absolute-coefficient polynomial at
// the right endpoint bounds |q''| everywhere in the interval. Positive bound
// arithmetic rounds upward; underflow is refused rather than reporting zero.
fn polynomial_error(
    evaluation: &FissionEnergyRelease,
    component: FissionEnergyComponent,
    a: f64,
    b: f64,
) -> Result<f64, EndfError> {
    let coefficients = evaluation.coefficients(component);
    if coefficients.len() < 3 {
        return Ok(0.0);
    }
    let high_ev = b / EV_J;
    let width_ev = (b - a) / EV_J;
    if !nonnegative(high_ev) || !nonnegative(width_ev) {
        return Err(FissionError::Arithmetic.into());
    }
    let mut derivative = 0.0;
    for k in (2..coefficients.len()).rev() {
        let coefficient =
            upper_product(coefficients[k].value.abs(), (k as f64) * ((k - 1) as f64))?;
        derivative = upper_sum(upper_product(derivative, high_ev)?, coefficient)?;
    }
    if derivative == 0.0 {
        return Ok(0.0);
    }
    let error_ev = upper_product(upper_product(derivative, width_ev * 0.5)?, width_ev * 0.25)?;
    upper_product(error_ev, EV_J)
}

#[allow(clippy::too_many_arguments)]
fn sher_beck_curve(
    evaluation: &FissionEnergyRelease,
    component: FissionEnergyComponent,
    domain: [f64; 2],
    neutron_yield: &FissionNeutronYield,
    reference_nu: f64,
    work: &mut Work,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<TabulatedYield, EndfError> {
    let base = evaluation.reference_value_j(component) / EV_J;
    let (slope, uses_nu): (f64, bool) = match component {
        FissionEnergyComponent::PromptNeutrons => (1.307, true),
        FissionEnergyComponent::Recoverable => (1.157, true),
        FissionEnergyComponent::Total => (1.057, true),
        FissionEnergyComponent::DelayedPhotons | FissionEnergyComponent::Betas => (-0.075, false),
        FissionEnergyComponent::Neutrinos => (-0.100, false),
        _ => (0.0, false),
    };
    let value = |energy_j: f64, nu: f64| -> Result<f64, EndfError> {
        let correction = if uses_nu {
            -8.07e6 * (nu - reference_nu)
        } else {
            0.0
        };
        let energy = slope.mul_add(energy_j / EV_J, base) + correction;
        if nonnegative(energy) {
            Ok(energy)
        } else {
            Err(invalid(
                "Sher-Beck energy is negative or unrepresentable on the requested domain",
            ))
        }
    };
    if !uses_nu {
        work.evaluate(checkpoint)?;
        let left = value(domain[0], 0.0)?;
        work.evaluate(checkpoint)?;
        let right = value(domain[1], 0.0)?;
        work.points(2)?;
        return make_curve(
            evaluation,
            component,
            domain.to_vec(),
            vec![left, right],
            checkpoint,
        );
    }
    let curve = neutron_yield.curve();
    let mut energies = Vec::new();
    let mut values = Vec::new();
    let mut region = 0;
    for (index, pair) in curve.energies_j().windows(2).enumerate() {
        checkpoint()?;
        if pair[0] == pair[1] || pair[1] <= domain[0] || pair[0] >= domain[1] {
            continue;
        }
        while curve.interpolation_regions()[region].end_point <= index + 1 {
            region += 1;
        }
        let a = pair[0].max(domain[0]);
        let b = pair[1].min(domain[1]);
        let mut endpoints = [0.0; 2];
        for (output, energy) in endpoints.iter_mut().zip([a, b]) {
            work.evaluate(checkpoint)?;
            let nu = match curve.interpolation_regions()[region].interpolation {
                EndfInterpolation::Histogram => curve.values()[index],
                EndfInterpolation::Linear => {
                    let width = pair[1] - pair[0];
                    convex(
                        curve.values()[index],
                        curve.values()[index + 1],
                        (pair[1] - energy) / width,
                        (energy - pair[0]) / width,
                    )?
                }
            };
            *output = value(energy, nu)?;
        }
        if energies.last().copied() != Some(a) || values.last().copied() != Some(endpoints[0]) {
            work.points(1)?;
            energies.push(a);
            values.push(endpoints[0]);
        }
        work.points(1)?;
        energies.push(b);
        values.push(endpoints[1]);
    }
    work.evaluate(checkpoint)?;
    let endpoint = value(domain[1], curve.value_at_energy(domain[1])?)?;
    if values.last().copied() != Some(endpoint) {
        work.points(1)?;
        energies.push(domain[1]);
        values.push(endpoint);
    }
    make_curve(evaluation, component, energies, values, checkpoint)
}

fn upper_product(a: f64, b: f64) -> Result<f64, EndfError> {
    if a == 0.0 || b == 0.0 {
        return Ok(0.0);
    }
    let value = a * b;
    if !value.is_finite() || value <= 0.0 || !value.next_up().is_finite() {
        Err(FissionError::Arithmetic.into())
    } else {
        Ok(value.next_up())
    }
}

fn upper_sum(a: f64, b: f64) -> Result<f64, EndfError> {
    if a == 0.0 {
        return Ok(b);
    }
    if b == 0.0 {
        return Ok(a);
    }
    let value = a + b;
    if !value.is_finite() || !value.next_up().is_finite() {
        Err(FissionError::Arithmetic.into())
    } else {
        Ok(value.next_up())
    }
}

fn to_joules(value_ev: f64) -> Result<f64, EndfError> {
    let joules = value_ev * EV_J;
    if !nonnegative(joules) || (value_ev > 0.0 && joules == 0.0) {
        Err(FissionError::Arithmetic.into())
    } else {
        Ok(joules)
    }
}

fn invalid(message: &'static str) -> EndfError {
    FissionError::InvalidInput(message).into()
}
