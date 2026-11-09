//! Source-bound, finite-capacity immobile storage on the original aqueous fluxes.
//!
//! The original `AqueousLangmuir` owns q(c), its vacant fraction and tangent.
//! An independent `LocalLangmuirTransfer` optionally supplies the local LDF
//! law dq/dt=k(q_eq(c)-q). Backward Euler eliminates loading analytically;
//! only the concentration storage becomes nonlinear. This is a first-order
//! passive analytical-element model, not aqueous speciation or site competition.

use super::sorption::{SorbingDiluteSoluteStep, finish_sorption};
use super::transport::{
    self, DiluteSoluteControls, DiluteSoluteError, DiluteSoluteStep, ImmobileStorage, Result,
    finite, integrated_flux, nonnegative, poll, positive,
};
use crate::aqueous_adsorption::{
    AdsorbentLoading, AqueousAdsorptionError, AqueousLangmuir, CONCENTRATION_AXIS,
    TEMPERATURE_AXIS,
    flow::{AdsorptionFlowError, LOADING_AXIS, LocalLangmuirTransfer, VOLUME_AXIS},
};
use fs_la::tridiagonal::TridiagonalWorkspace;
use fs_matdb::QueryPoint;
use fs_qty::{AmountConcentration, QuantitySpec, Temperature, Volume};

pub(super) mod populations;

/// Original equilibrium law or its independently sourced LOCAL kinetic law.
/// A finite-batch clock/volume correlation cannot substitute for the local law.
#[derive(Clone, Copy, Debug)]
pub enum LangmuirSoluteLaw<'a> {
    Equilibrium(&'a AqueousLangmuir),
    Kinetic(&'a LocalLangmuirTransfer),
}
impl<'a> LangmuirSoluteLaw<'a> {
    fn equilibrium(self) -> &'a AqueousLangmuir {
        match self {
            Self::Equilibrium(m) => m,
            Self::Kinetic(m) => m.equilibrium(),
        }
    }
    fn rate(self) -> Option<f64> {
        match self {
            Self::Equilibrium(_) => None,
            Self::Kinetic(m) => Some(m.rate().value()),
        }
    }
    fn limits(self) -> [f64; 2] {
        let (mut low, mut high) = self
            .equilibrium()
            .source()
            .validity
            .bound(CONCENTRATION_AXIS)
            .expect("admitted original concentration support");
        if let Self::Kinetic(m) = self {
            let (a, b) = m
                .source()
                .validity
                .bound(CONCENTRATION_AXIS)
                .expect("admitted local concentration support");
            low = low.max(a);
            high = high.min(b);
        }
        [low, high]
    }
}

/// One independent sorbent population for a cell/constituent slot. `None` in
/// the input array denotes an inert slot with zero occupied and vacant stocks.
#[derive(Clone, Copy, Debug)]
pub struct LangmuirSoluteSorption<'a> {
    pub law: LangmuirSoluteLaw<'a>,
    /// Original retained dry adsorbent mass per bulk transverse area [kg/m²].
    pub dry_mass_kg_m2: f64,
    /// Prescribed physical temperature; no sorption heat or thermal feedback.
    pub temperature: Temperature,
    /// Additional original conditions only. Actual concentration, loading,
    /// temperature and cell solution volume cannot be overridden here.
    pub conditions: &'a QueryPoint,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LangmuirSoluteControls<const N: usize> {
    pub transport: DiluteSoluteControls<N>,
    /// Newton iterations for each constituent, in 1..=256.
    pub maximum_iterations: usize,
    /// Box-preserving residual backtracks per iteration, in 1..=64.
    pub maximum_line_searches: usize,
    /// Cumulative native tridiagonal solves for all constituents, in 1..=8192.
    pub maximum_linear_solves: usize,
    /// Every actual local source point, including rejected trials, in 1..=10^7.
    pub maximum_source_evaluations: usize,
}

/// Original source errors retain their typed cause and exact cell/species.
#[derive(Clone, Debug, PartialEq)]
pub enum LangmuirSoluteError {
    Transport(DiluteSoluteError),
    EquilibriumSource {
        cell: usize,
        species: usize,
        error: AqueousAdsorptionError,
    },
    KineticSource {
        cell: usize,
        species: usize,
        error: AdsorptionFlowError,
    },
    SourceBudget {
        evaluations: usize,
        maximum_evaluations: usize,
    },
}
impl From<DiluteSoluteError> for LangmuirSoluteError {
    fn from(value: DiluteSoluteError) -> Self {
        Self::Transport(value)
    }
}
impl core::fmt::Display for LangmuirSoluteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "porous Langmuir solute transport: {self:?}")
    }
}
impl std::error::Error for LangmuirSoluteError {}

/// Accepted two-phase transport plus independently retained unoccupied sites.
#[derive(Clone, Debug, PartialEq)]
pub struct LangmuirSoluteStep<const N: usize> {
    pub sorption: SorbingDiluteSoluteStep<N>,
    /// Vacant sites [mol/m²] are not solute and never enter the metal balance.
    /// Supply this field, alongside the occupied stock, to the next interval.
    pub vacant_amount_mol_m2: Vec<[f64; N]>,
    pub maximum_site_balance_defect_mol_m2: [f64; N],
    pub source_evaluations: usize,
}

pub(super) struct LangmuirStorage<'a, const N: usize> {
    bindings: &'a [[Option<LangmuirSoluteSorption<'a>>; N]],
    original: &'a [[f64; N]],
    original_vacant: &'a [[f64; N]],
    area: f64,
    duration: f64,
    controls: LangmuirSoluteControls<N>,
    initial_concentration: Vec<[f64; N]>,
    final_water: Vec<f64>,
    bound: Vec<[f64; N]>,
    vacant: Vec<[f64; N]>,
    exchange: Vec<[f64; N]>,
    maximum_defect: [f64; N],
    maximum_site_defect: [f64; N],
    source_evaluations: usize,
    linear_solves: usize,
    source_error: Option<LangmuirSoluteError>,
}

impl<'a, const N: usize> LangmuirStorage<'a, N> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        n: usize,
        duration: f64,
        original: &'a [[f64; N]],
        original_vacant: &'a [[f64; N]],
        bindings: &'a [[Option<LangmuirSoluteSorption<'a>>; N]],
        area: f64,
        controls: LangmuirSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<Self, LangmuirSoluteError> {
        poll(cancelled)?;
        let (_, base) = transport::admit_extent(n, controls.transport)?;
        let bytes = base + Self::extra_bytes(n);
        if bytes > controls.transport.maximum_storage_bytes {
            return Err(DiluteSoluteError::Workspace {
                required_bytes: bytes,
                maximum_bytes: controls.transport.maximum_storage_bytes,
            }
            .into());
        }
        if original.len() != n
            || original_vacant.len() != n
            || bindings.len() != n
            || !positive(area)
            || !positive(duration)
            || !(1..=256).contains(&controls.maximum_iterations)
            || !(1..=64).contains(&controls.maximum_line_searches)
            || !(1..=8192).contains(&controls.maximum_linear_solves)
            || !(1..=10_000_000).contains(&controls.maximum_source_evaluations)
        {
            return Err(DiluteSoluteError::Invalid("matching Langmuir extents, positive physical area/time and bounded nonlinear work required").into());
        }
        let mut elements = [None; N];
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            for k in 0..N {
                let (old, vacant) = (original[i][k], original_vacant[i][k]);
                if !nonnegative(old) || !nonnegative(vacant) {
                    return Err(DiluteSoluteError::Invalid(
                        "finite nonnegative occupied and vacant stocks required",
                    )
                    .into());
                }
                if let Some(b) = bindings[i][k] {
                    if !positive(b.dry_mass_kg_m2)
                        || !positive(b.temperature.value())
                        || b.law.limits()[0] > b.law.limits()[1]
                        || [
                            CONCENTRATION_AXIS,
                            TEMPERATURE_AXIS,
                            LOADING_AXIS,
                            VOLUME_AXIS,
                        ]
                        .iter()
                        .any(|axis| b.conditions.axes().contains_key(*axis))
                    {
                        return Err(DiluteSoluteError::Invalid("positive dry mass/temperature, intersecting sources and unmodified actual coordinates required").into());
                    }
                    let element = b.law.equilibrium().element();
                    if elements[k].is_some_and(|e| e != element) {
                        return Err(DiluteSoluteError::Invalid("every cell of a transported constituent must name the same analytical element").into());
                    }
                    elements[k] = Some(element);
                    let capacity =
                        integrated_flux(b.dry_mass_kg_m2, b.law.equilibrium().capacity().value())?;
                    if old > capacity
                        || vacant > capacity
                        || site_balance(old, vacant, capacity)?.abs()
                            > controls.transport.absolute_equation_tolerance_mol_m2[k]
                    {
                        return Err(DiluteSoluteError::Invalid(
                            "occupied plus vacant stocks must equal original finite site capacity",
                        )
                        .into());
                    }
                } else if old != 0.0 || vacant != 0.0 {
                    return Err(DiluteSoluteError::Invalid(
                        "an inert constituent slot cannot retain occupied or vacant sites",
                    )
                    .into());
                }
            }
        }
        for k in 0..N {
            if elements[k].is_some() && elements[..k].contains(&elements[k]) {
                return Err(DiluteSoluteError::Invalid(
                    "duplicate analytical element in independent constituent slots",
                )
                .into());
            }
        }
        Ok(Self {
            bindings,
            original,
            original_vacant,
            area,
            duration,
            controls,
            initial_concentration: vec![[0.0; N]; n],
            final_water: vec![0.0; n],
            bound: vec![[0.0; N]; n],
            vacant: vec![[0.0; N]; n],
            exchange: vec![[0.0; N]; n],
            maximum_defect: [0.0; N],
            maximum_site_defect: [0.0; N],
            source_evaluations: 0,
            linear_solves: 0,
            source_error: None,
        })
    }

    fn extra_bytes(n: usize) -> usize {
        // Four retained cell/species fields, final water, one species guess and
        // bounds, eight Newton arrays, and overlapping final phase reports.
        // Original source-owner query/receipt allocations are outside this
        // numerical payload allowance, as are the borrowed cards and inputs.
        (4 * n * N + 12 * n + 20 * N) * core::mem::size_of::<f64>() + 1024
    }

    fn refuse(&mut self, error: LangmuirSoluteError) -> DiluteSoluteError {
        self.source_error = Some(error);
        DiluteSoluteError::Invalid("original Langmuir source or evaluation budget refused")
    }

    pub(super) fn recover<T>(
        &mut self,
        value: Result<T>,
    ) -> core::result::Result<T, LangmuirSoluteError> {
        value.map_err(|error| {
            self.source_error
                .take()
                .unwrap_or(LangmuirSoluteError::Transport(error))
        })
    }

    fn tick(&mut self) -> Result<()> {
        if self.source_evaluations >= self.controls.maximum_source_evaluations {
            return Err(self.refuse(LangmuirSoluteError::SourceBudget {
                evaluations: self.source_evaluations,
                maximum_evaluations: self.controls.maximum_source_evaluations,
            }));
        }
        self.source_evaluations += 1;
        Ok(())
    }

    /// Complete actual endpoint state is checked by the original owner, after
    /// the same owner's constitutive arithmetic determines the implicit stock.
    fn admit_source(&mut self, i: usize, k: usize, c: f64, water: f64, bound: f64) -> Result<()> {
        let Some(b) = self.bindings[i][k] else {
            return Ok(());
        };
        admit_site_source(b, i, k, c, water, self.area, bound).map_err(|error| self.refuse(error))
    }

    fn response(&mut self, i: usize, k: usize, c: f64) -> Result<(f64, f64)> {
        let Some(b) = self.bindings[i][k] else {
            return Ok((0.0, 0.0));
        };
        self.tick()?;
        let response = site_response(
            b,
            i,
            k,
            self.original[i][k],
            self.original_vacant[i][k],
            c,
            self.final_water[i],
            self.area,
            self.duration,
            self.controls.transport,
        )
        .map_err(|error| match error {
            LangmuirSoluteError::Transport(error) => error,
            error => self.refuse(error),
        })?;
        self.bound[i][k] = response.bound;
        self.vacant[i][k] = response.vacant;
        self.exchange[i][k] = response.exchange;
        self.maximum_defect[k] = self.maximum_defect[k].max(response.equation_defect.abs());
        self.maximum_site_defect[k] = self.maximum_site_defect[k].max(response.site_defect.abs());
        Ok((response.exchange, response.tangent))
    }

    pub(super) fn finish(
        self,
        initial: &[[f64; N]],
        step: DiluteSoluteStep<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<LangmuirSoluteStep<N>, LangmuirSoluteError> {
        let sorption = finish_sorption(
            self.original,
            self.bound,
            self.exchange,
            self.maximum_defect,
            initial,
            step,
            cancelled,
        )?;
        Ok(LangmuirSoluteStep {
            sorption,
            vacant_amount_mol_m2: self.vacant,
            maximum_site_balance_defect_mol_m2: self.maximum_site_defect,
            source_evaluations: self.source_evaluations,
        })
    }
}

impl<const N: usize> ImmobileStorage<N> for LangmuirStorage<'_, N> {
    fn additional_bytes(&self) -> usize {
        Self::extra_bytes(self.original.len())
    }

    fn admit(
        &mut self,
        mobile: &[[f64; N]],
        before: &[f64],
        after: &[f64],
        controls: DiluteSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        for i in 0..mobile.len() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            self.final_water[i] = after[i];
            for k in 0..N {
                let c = finite(mobile[i][k] / before[i])?;
                if mobile[i][k] > 0.0 && c == 0.0 {
                    return Err(DiluteSoluteError::Invalid(
                        "positive initial concentration is unrepresentable",
                    ));
                }
                self.initial_concentration[i][k] = c;
                if let Some(b) = self.bindings[i][k] {
                    self.tick()?;
                    self.admit_source(i, k, c, before[i], self.original[i][k])?;
                    if matches!(b.law, LangmuirSoluteLaw::Equilibrium(_)) {
                        let (q, vacant_q, _) =
                            b.law.equilibrium().loading_components(c).map_err(|error| {
                                self.refuse(LangmuirSoluteError::EquilibriumSource {
                                    cell: i,
                                    species: k,
                                    error,
                                })
                            })?;
                        for residual in [
                            self.original[i][k] - integrated_flux(b.dry_mass_kg_m2, q)?,
                            self.original_vacant[i][k]
                                - integrated_flux(b.dry_mass_kg_m2, vacant_q)?,
                        ] {
                            if residual.abs() > controls.absolute_equation_tolerance_mol_m2[k] {
                                return Err(DiluteSoluteError::Sorption {
                                    cell: i,
                                    species: k,
                                    residual_mol_m2: residual,
                                    allowance_mol_m2: controls.absolute_equation_tolerance_mol_m2
                                        [k],
                                });
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn solve_species(
        &mut self,
        k: usize,
        lower: &[f64],
        diagonal: &[f64],
        upper: &[f64],
        rhs: &mut [f64],
        solver: &mut TridiagonalWorkspace,
        controls: DiluteSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize> {
        if self
            .bindings
            .iter()
            .all(|row| row[k].is_none_or(|b| b.law.rate() == Some(0.0)))
        {
            if self.linear_solves == self.controls.maximum_linear_solves {
                return Err(DiluteSoluteError::NonlinearWork {
                    required_linear_solves: self.linear_solves + 1,
                    maximum_linear_solves: self.controls.maximum_linear_solves,
                });
            }
            solver.factor(lower, diagonal, upper, &mut *cancelled)?;
            solver.solve_factored_in_place(rhs, &mut *cancelled)?;
            self.linear_solves += 1;
            return Ok(1);
        }
        let initial: Vec<_> = self.initial_concentration.iter().map(|c| c[k]).collect();
        let limits: Vec<_> = self
            .bindings
            .iter()
            .map(|row| row[k].map_or([0.0, f64::MAX], |b| b.law.limits()))
            .collect();
        let tolerance = controls.absolute_equation_tolerance_mol_m2[k]
            .min(controls.amount_balance_tolerance_mol_m2[k] / rhs.len() as f64)
            * 0.25;
        let ctl = transport::nonlinear::Controls {
            tolerance_mol_m2: tolerance,
            maximum_iterations: self.controls.maximum_iterations,
            maximum_line_searches: self.controls.maximum_line_searches,
            remaining_linear_solves: self.controls.maximum_linear_solves - self.linear_solves,
        };
        let count = transport::nonlinear::solve(
            k,
            lower,
            diagonal,
            upper,
            rhs,
            &initial,
            |i| limits[i],
            solver,
            ctl,
            |i, c| self.response(i, k, c),
            cancelled,
        )?;
        self.linear_solves += count;
        Ok(count)
    }

    fn finish_cell(
        &mut self,
        i: usize,
        k: usize,
        c: f64,
        _: DiluteSoluteControls<N>,
    ) -> Result<()> {
        if i == 0 {
            // Trial defects are private; final reports describe accepted cells.
            self.maximum_defect[k] = 0.0;
            self.maximum_site_defect[k] = 0.0;
        }
        self.response(i, k, c).map(|_| ())
    }
    fn amounts(&self, i: usize, k: usize) -> [f64; 3] {
        [self.original[i][k], self.bound[i][k], self.exchange[i][k]]
    }
}

/// One physical site's response, shared by the original one-population path
/// and the sparse multiple-population adapter. Source admission and arithmetic
/// stay here; neither spatial adapter invents its own isotherm or kinetic law.
struct SiteResponse {
    bound: f64,
    vacant: f64,
    exchange: f64,
    tangent: f64,
    equation_defect: f64,
    site_defect: f64,
}

#[allow(clippy::too_many_arguments)]
fn admit_site_source(
    b: LangmuirSoluteSorption<'_>,
    i: usize,
    k: usize,
    c: f64,
    water: f64,
    area: f64,
    bound: f64,
) -> core::result::Result<(), LangmuirSoluteError> {
    let loading = finite(bound / b.dry_mass_kg_m2)?;
    if bound > 0.0 && loading == 0.0 {
        return Err(DiluteSoluteError::Invalid("positive loading is unrepresentable").into());
    }
    let volume = Volume::new(integrated_flux(area, water)?);
    match b.law {
        LangmuirSoluteLaw::Equilibrium(model) => {
            let point = b
                .conditions
                .clone()
                .with_quantity(
                    LOADING_AXIS,
                    QuantitySpec::dimensional(AdsorbentLoading::DIMS),
                    loading,
                )
                .and_then(|p| {
                    p.with_quantity(
                        VOLUME_AXIS,
                        QuantitySpec::dimensional(Volume::DIMS),
                        volume.value(),
                    )
                })
                .map_err(|e| LangmuirSoluteError::EquilibriumSource {
                    cell: i,
                    species: k,
                    error: e.into(),
                })?;
            model
                .admit_loading_point(AmountConcentration::new(c), b.temperature, &point)
                .map_err(|error| LangmuirSoluteError::EquilibriumSource {
                    cell: i,
                    species: k,
                    error,
                })?;
        }
        LangmuirSoluteLaw::Kinetic(model) => {
            model
                .admit_loading_point(
                    AmountConcentration::new(c),
                    AdsorbentLoading::new(loading),
                    volume,
                    b.temperature,
                    b.conditions,
                )
                .map_err(|error| LangmuirSoluteError::KineticSource {
                    cell: i,
                    species: k,
                    error,
                })?;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn site_response<const N: usize>(
    b: LangmuirSoluteSorption<'_>,
    i: usize,
    k: usize,
    original: f64,
    original_vacant: f64,
    c: f64,
    water: f64,
    area: f64,
    duration: f64,
    controls: DiluteSoluteControls<N>,
) -> core::result::Result<SiteResponse, LangmuirSoluteError> {
    let (retained, weight) = weights(b.law.rate(), duration)?;
    if weight == 0.0 {
        admit_site_source(b, i, k, c, water, area, original)?;
        return Ok(SiteResponse {
            bound: original,
            vacant: original_vacant,
            exchange: 0.0,
            tangent: 0.0,
            equation_defect: 0.0,
            site_defect: 0.0,
        });
    }
    let (q, vacant_q, tangent_q) = b.law.equilibrium().loading_components(c).map_err(|error| {
        LangmuirSoluteError::EquilibriumSource {
            cell: i,
            species: k,
            error,
        }
    })?;
    let target = integrated_flux(b.dry_mass_kg_m2, q)?;
    let vacant_target = integrated_flux(b.dry_mass_kg_m2, vacant_q)?;
    let endpoint = interpolate(original, target, retained, weight)?;
    let vacant = interpolate(original_vacant, vacant_target, retained, weight)?;
    admit_site_source(b, i, k, c, water, area, endpoint)?;
    let actual = finite(endpoint - original)?;
    let target_difference = if original > original_vacant && target > vacant_target {
        finite(original_vacant - vacant_target)?
    } else {
        finite(target - original)?
    };
    let predicted = integrated_flux(weight, target_difference)?;
    let residual = finite(actual - predicted)?;
    let allowance = finite(
        controls.absolute_equation_tolerance_mol_m2[k]
            + controls.relative_equation_tolerance * actual.abs().max(predicted.abs()),
    )?;
    let capacity = integrated_flux(b.dry_mass_kg_m2, b.law.equilibrium().capacity().value())?;
    let site_defect = site_balance(endpoint, vacant, capacity)?;
    if endpoint > capacity
        || vacant > capacity
        || residual.abs() > allowance
        || site_defect.abs() > controls.absolute_equation_tolerance_mol_m2[k]
    {
        return Err(DiluteSoluteError::Sorption {
            cell: i,
            species: k,
            residual_mol_m2: if residual.abs() > allowance {
                residual
            } else {
                site_defect
            },
            allowance_mol_m2: allowance,
        }
        .into());
    }
    Ok(SiteResponse {
        bound: endpoint,
        vacant,
        exchange: actual,
        tangent: integrated_flux(weight, integrated_flux(b.dry_mass_kg_m2, tangent_q)?)?,
        equation_defect: residual,
        site_defect,
    })
}

fn weights(rate: Option<f64>, dt: f64) -> Result<(f64, f64)> {
    let Some(k) = rate else {
        return Ok((0.0, 1.0));
    };
    let x = k * dt;
    if !x.is_finite() {
        return Err(DiluteSoluteError::Invalid(
            "finite Langmuir rate-times-duration required",
        ));
    }
    let (retained, weight) = if x <= 1.0 {
        (1.0 / (1.0 + x), x / (1.0 + x))
    } else {
        let inverse = (1.0 / k) / dt;
        (inverse / (1.0 + inverse), 1.0 / (1.0 + inverse))
    };
    if (k > 0.0 && weight == 0.0) || !retained.is_finite() || !weight.is_finite() {
        return Err(DiluteSoluteError::Invalid(
            "Langmuir relaxation fractions are unrepresentable",
        ));
    }
    Ok((retained, weight))
}

fn interpolate(old: f64, target: f64, retained: f64, weight: f64) -> Result<f64> {
    if old == target {
        return Ok(old);
    }
    let (base, other, fraction) = if weight <= retained {
        (old, target, weight)
    } else {
        (target, old, retained)
    };
    let result = finite(base + integrated_flux(fraction, finite(other - base)?)?)?;
    if result < old.min(target) || result > old.max(target) {
        return Err(DiluteSoluteError::Invalid(
            "implicit Langmuir stock lost its convex range",
        ));
    }
    Ok(result)
}

fn site_balance(occupied: f64, vacant: f64, capacity: f64) -> Result<f64> {
    compensated_sum([occupied, vacant, -capacity])
}

fn compensated_sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = finite(total + value)?;
        correction += if total.abs() >= value.abs() {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    finite(total + correction)
}

#[cfg(test)]
mod tests;
