//! Finite aqueous adsorption/desorption and diffusion on original surface sites.
//!
//! The original `AqueousSolution` supplies the finite bulk inventory and its
//! explicit convex `VirialSolutionModel`, including nonideal co-solute effects.
//! A zero interaction matrix is the ideal limit; a datum-free ideal inventory
//! is deliberately not assigned an invented chemical reference. All surface
//! faces represent one explicitly selected, nonionic original bulk species.
//!
//! For bulk activity `a`, surface occupancy `theta` and the authored standard
//! potential difference `delta_mu = mu_surface^0 - mu_bulk^0`, detailed balance
//! requires `theta/(1-theta) = a*exp(-delta_mu/RT)`. The original scalar site
//! owner therefore receives adsorption rate `kd*a*exp(-delta_mu/RT)` and
//! desorption rate `kd`. A common offset across the surface also leaves the
//! existing occupancy-driven partition diffusion thermodynamically consistent.
//! No adsorption affinity or rate is inferred from surface tension.
//!
//! All sites exchange simultaneously through a monotone scalar finite-bath
//! balance. Their original backward-Euler maps preserve finite packing; the
//! same bulk amount is debited once. Half surface-diffusion sweeps surround
//! each exchange step. Backward Euler makes the complete composition FIRST
//! order in time, despite the symmetric diffusion arrangement. Refine the
//! explicit maximum step; no temporal-error or global accuracy certificate is
//! claimed. Areas, temperature, solvent volume, topology and coefficients stay
//! fixed. There is no bulk spatial transport, micellization, ionic adsorption,
//! adsorption heat or mechanical motion in this adapter.
//!
//! Original amount/energy laws, declared concentration support, work bounds and
//! cancellation are checked before either original inventory is published.
//! All supplied coefficients retain caller-owned material applicability; raw
//! numbers establish no source authority. Changing the bulk reference by a
//! factor `s` requires changing `delta_mu` by `-RT*ln(s)` for the same physics.

use super::diffusion::{
    SurfaceDiffusionControls, SurfaceDiffusionError, SurfacePartitionDiffusionEdge,
    diffuse_partitioned_surface_surfactant,
};
use super::{Sum, SurfactantError, SurfactantFace, SurfactantLaw, det};
use fs_material::adsorption::competitive::sites::scalar::{ScalarSiteError, backward_euler_site};
use fs_material::osmosis::{AqueousSolution, MAX_OSMOTIC_SOLUTES, OsmosisError};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Original species layout and explicit kinetic/chemical coupling.
#[derive(Clone, Copy, Debug)]
pub struct SolubleSurfaceKinetics<'a> {
    /// Index in the original solution's caller-owned independent-species order.
    pub bulk_species: usize,
    /// Surface-minus-bulk standard chemical potential on the original bulk
    /// model's activity reference [J/mol], common to every surface face.
    pub standard_potential_offset_j_mol: f64,
    /// One held desorption rate per face [1/s]. Zero seals that face to bulk
    /// exchange, while its original lateral surface diffusion remains active.
    pub desorption_rates_s_inverse: &'a [f64],
}

/// Explicit whole-call work, physical-balance and time-resolution limits.
#[derive(Clone, Copy, Debug)]
pub struct SolubleSurfaceControls {
    pub maximum_step_s: f64,
    pub maximum_substeps: usize,
    /// Every scalar-balance evaluation, including bracket endpoints.
    pub maximum_root_evaluations: usize,
    /// Pair applications in both surface-diffusion halves of every substep.
    pub maximum_surface_pair_exchanges: usize,
    pub maximum_faces: usize,
    pub maximum_edges: usize,
    /// Absolute residual of the original finite-bath amount equation [mol].
    pub root_tolerance_mol: f64,
    /// Original-to-final combined bulk plus surface amount allowance [mol].
    pub amount_tolerance_mol: f64,
    /// Original-to-final and per-stage free-energy increase allowance [J].
    /// Also bounds discrepancy between endpoint and accumulated stage changes.
    pub energy_tolerance_j: f64,
}

/// Physical transfers and the independently reconstructed closed-system account.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SolubleSurfaceReport {
    pub duration_s: f64,
    pub substeps: usize,
    pub root_evaluations: usize,
    pub surface_pair_exchanges: usize,
    pub bulk_amount_before_mol: f64,
    pub bulk_amount_after_mol: f64,
    pub surface_amount_before_mol: f64,
    pub surface_amount_after_mol: f64,
    /// Integrated original BE adsorption and desorption streams, not net uptake.
    pub gross_adsorption_mol: f64,
    pub gross_desorption_mol: f64,
    pub maximum_root_residual_mol: f64,
    pub maximum_site_kinetic_residual_mol: f64,
    pub amount_defect_mol: f64,
    /// Original bulk mixing plus original surface energy plus offset*N_surface.
    pub free_energy_change_j: f64,
    pub exchange_free_energy_dissipation_j: f64,
    pub surface_diffusion_free_energy_dissipation_j: f64,
    pub energy_account_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SolubleSurfaceError {
    Invalid(&'static str),
    Material(OsmosisError),
    Surface { face: usize, error: SurfactantError },
    Site { face: usize, error: ScalarSiteError },
    Diffusion(SurfaceDiffusionError),
    Allocation,
    Unrepresentable,
    WorkBudget,
    RootBudget { completed_evaluations: usize },
    BulkSupport { residual_mol: f64 },
    RootAccuracy { residual_mol: f64 },
    AmountBalance { defect_mol: f64 },
    FreeEnergy { change_j: f64 },
    EnergyAccount { residual_j: f64 },
    Cancelled,
}
impl core::fmt::Display for SolubleSurfaceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "soluble surface transport: {self:?}")
    }
}
impl std::error::Error for SolubleSurfaceError {}
impl From<OsmosisError> for SolubleSurfaceError {
    fn from(error: OsmosisError) -> Self {
        Self::Material(error)
    }
}
impl From<SurfaceDiffusionError> for SolubleSurfaceError {
    fn from(error: SurfaceDiffusionError) -> Self {
        Self::Diffusion(error)
    }
}
type Result<T> = core::result::Result<T, SolubleSurfaceError>;

#[derive(Clone, Copy)]
struct Site {
    capacity: f64,
    desorption: f64,
    log_adsorption_prefactor: f64,
}

/// Advance one finite bath and the same surface face amounts atomically.
///
/// The source model, original species order and all co-solute inventories are
/// retained. Every face must use its original Langmuir EOS at the exact bulk
/// temperature. Negative tension, support exhaustion, unrepresentable rates,
/// exhausted work, physical-account failure and cancellation publish nothing.
#[allow(clippy::too_many_arguments)]
pub fn advance_soluble_surface(
    bulk: &mut AqueousSolution,
    laws: &[SurfactantLaw],
    areas_m2: &[f64],
    amounts_mol: &mut [f64],
    kinetics: SolubleSurfaceKinetics<'_>,
    edges: &[SurfacePartitionDiffusionEdge],
    duration_s: f64,
    controls: SolubleSurfaceControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SolubleSurfaceReport> {
    poll(&mut cancelled)?;
    let species = kinetics.bulk_species;
    let count = bulk.solute_moles().len();
    bulk.nonideal_model().ok_or(SolubleSurfaceError::Invalid(
        "bulk exchange requires the original explicit activity/free-energy model",
    ))?;
    if species >= count
        || count > MAX_OSMOTIC_SOLUTES
        || areas_m2.is_empty()
        || laws.len() != areas_m2.len()
        || amounts_mol.len() != areas_m2.len()
        || kinetics.desorption_rates_s_inverse.len() != areas_m2.len()
        || !kinetics.standard_potential_offset_j_mol.is_finite()
        || !duration_s.is_finite()
        || duration_s < 0.0
        || !controls.maximum_step_s.is_finite()
        || controls.maximum_step_s <= 0.0
        || controls.maximum_substeps == 0
        || controls.maximum_root_evaluations == 0
        || !controls.root_tolerance_mol.is_finite()
        || controls.root_tolerance_mol <= 0.0
        || !nonnegative(controls.amount_tolerance_mol)
        || !nonnegative(controls.energy_tolerance_j)
    {
        return Err(SolubleSurfaceError::Invalid(
            "matching original layouts and finite physical/time/work controls required",
        ));
    }
    if areas_m2.len() > controls.maximum_faces || edges.len() > controls.maximum_edges {
        return Err(SolubleSurfaceError::WorkBudget);
    }
    let rt = finite(R * bulk.temperature_k())?;
    let offset = kinetics.standard_potential_offset_j_mol;
    let mut sites = Vec::new();
    sites
        .try_reserve_exact(areas_m2.len())
        .map_err(|_| SolubleSurfaceError::Allocation)?;
    for (face, ((&law, &area), &rate)) in laws
        .iter()
        .zip(areas_m2)
        .zip(kinetics.desorption_rates_s_inverse)
        .enumerate()
    {
        poll_at(&mut cancelled, face)?;
        let SurfactantLaw::Langmuir {
            thermal_energy_j_mol,
            saturation_concentration_mol_m2,
            ..
        } = law
        else {
            return Err(SolubleSurfaceError::Invalid(
                "original Langmuir sites required for soluble exchange",
            ));
        };
        if thermal_energy_j_mol != rt || !nonnegative(rate) {
            return Err(SolubleSurfaceError::Invalid(
                "one held temperature and finite nonnegative desorption rates required",
            ));
        }
        SurfactantFace {
            law,
            amount_mol: amounts_mol[face],
        }
        .at_area(area)
        .map_err(|error| SolubleSurfaceError::Surface { face, error })?;
        let capacity = finite(area * saturation_concentration_mol_m2)?;
        if capacity <= 0.0 {
            return Err(SolubleSurfaceError::Unrepresentable);
        }
        sites.push(Site {
            capacity,
            desorption: rate,
            log_adsorption_prefactor: if rate == 0.0 {
                0.0
            } else {
                finite(det::ln(rate) - offset / rt)?
            },
        });
    }
    let substeps = if duration_s == 0.0 {
        0
    } else {
        let required = (duration_s / controls.maximum_step_s).ceil().max(1.0);
        if !required.is_finite() || required > controls.maximum_substeps as f64 {
            return Err(SolubleSurfaceError::WorkBudget);
        }
        let n = required as usize;
        if n == 0 || n > controls.maximum_substeps {
            return Err(SolubleSurfaceError::WorkBudget);
        }
        n
    };
    let pair_work = if edges.is_empty() {
        0
    } else {
        edges
            .len()
            .checked_mul(2)
            .and_then(|x| x.checked_sub(1))
            .and_then(|x| x.checked_mul(2))
            .and_then(|x| x.checked_mul(substeps))
            .ok_or(SolubleSurfaceError::WorkBudget)?
    };
    if pair_work > controls.maximum_surface_pair_exchanges {
        return Err(SolubleSurfaceError::WorkBudget);
    }
    let step = if substeps == 0 {
        0.0
    } else {
        duration_s / substeps as f64
    };
    if substeps > 0 && (step <= 0.0 || (!edges.is_empty() && 0.5 * step == 0.0)) {
        return Err(SolubleSurfaceError::Unrepresentable);
    }
    let mut surface = copied(amounts_mol, &mut cancelled)?;
    let mut trial = copied(amounts_mol, &mut cancelled)?;
    let mut inventory = [0.0; MAX_OSMOTIC_SOLUTES];
    inventory[..count].copy_from_slice(bulk.solute_moles());
    let diffusion_controls = SurfaceDiffusionControls {
        maximum_step_s: controls.maximum_step_s,
        maximum_substeps: 1,
        maximum_faces: controls.maximum_faces,
        maximum_edges: controls.maximum_edges,
        amount_tolerance_mol: controls.amount_tolerance_mol,
        energy_increase_tolerance_j: controls.energy_tolerance_j,
    };
    // Admit the complete original surface graph even for a zero-duration call.
    let initial = diffuse_partitioned_surface_surfactant(
        laws,
        areas_m2,
        &mut surface,
        edges,
        0.0,
        diffusion_controls,
        &mut cancelled,
    )?;
    let mut report = SolubleSurfaceReport {
        duration_s,
        substeps,
        surface_pair_exchanges: pair_work,
        bulk_amount_before_mol: inventory[species],
        surface_amount_before_mol: initial.amount_before_mol,
        ..SolubleSurfaceReport::default()
    };
    let mut exchange_change = Sum::default();
    let mut diffusion_change = Sum::default();
    let mut adsorbed = Sum::default();
    let mut desorbed = Sum::default();
    for _ in 0..substeps {
        poll(&mut cancelled)?;
        for before_exchange in [true, false] {
            if !edges.is_empty() {
                let diffusion = diffuse_partitioned_surface_surfactant(
                    laws,
                    areas_m2,
                    &mut surface,
                    edges,
                    0.5 * step,
                    diffusion_controls,
                    &mut cancelled,
                )?;
                let change =
                    finite(diffusion.free_energy_change_j + offset * diffusion.amount_defect_mol)?;
                admit_energy(change, controls)?;
                diffusion_change.add(change);
            }
            if !before_exchange {
                continue;
            }
            let equation = Exchange {
                bulk,
                species,
                sites: &sites,
                original: &surface,
                bulk_amount: inventory[species],
                duration_s: step,
            };
            let exchange = equation.solve(
                &mut trial,
                controls,
                &mut report.root_evaluations,
                &mut cancelled,
            )?;
            let mut after = inventory;
            after[species] = exchange.bulk_amount;
            let change = energy_change(
                bulk,
                species,
                &inventory[..count],
                &after[..count],
                laws,
                areas_m2,
                areas_m2,
                &surface,
                &trial,
                offset,
                &mut cancelled,
            )?;
            admit_energy(change, controls)?;
            exchange_change.add(change);
            adsorbed.add(exchange.adsorption);
            desorbed.add(exchange.desorption);
            report.maximum_root_residual_mol = report
                .maximum_root_residual_mol
                .max(exchange.residual.abs());
            report.maximum_site_kinetic_residual_mol = report
                .maximum_site_kinetic_residual_mol
                .max(exchange.kinetic_residual);
            inventory = after;
            core::mem::swap(&mut surface, &mut trial);
        }
    }
    // Every source support check precedes the original owner's atomic replace.
    bulk.osmotic_pressure_at(bulk.solvent_volume_m3(), &inventory[..count])?;
    let mut defect = Sum::default();
    let mut surface_total = Sum::default();
    defect.add(inventory[species] - bulk.solute_moles()[species]);
    for (face, (&before, &after)) in amounts_mol.iter().zip(&surface).enumerate() {
        poll_at(&mut cancelled, face)?;
        defect.add(after - before);
        surface_total.add(after);
    }
    report.amount_defect_mol = finite(defect.value())?;
    if report.amount_defect_mol.abs() > controls.amount_tolerance_mol {
        return Err(SolubleSurfaceError::AmountBalance {
            defect_mol: report.amount_defect_mol,
        });
    }
    report.free_energy_change_j = energy_change(
        bulk,
        species,
        bulk.solute_moles(),
        &inventory[..count],
        laws,
        areas_m2,
        areas_m2,
        amounts_mol,
        &surface,
        offset,
        &mut cancelled,
    )?;
    admit_energy(report.free_energy_change_j, controls)?;
    report.exchange_free_energy_dissipation_j = finite(-exchange_change.value())?;
    report.surface_diffusion_free_energy_dissipation_j = finite(-diffusion_change.value())?;
    let mut account = Sum::default();
    account.add(report.free_energy_change_j);
    account.add(report.exchange_free_energy_dissipation_j);
    account.add(report.surface_diffusion_free_energy_dissipation_j);
    report.energy_account_residual_j = finite(account.value())?;
    if report.energy_account_residual_j.abs() > controls.energy_tolerance_j {
        return Err(SolubleSurfaceError::EnergyAccount {
            residual_j: report.energy_account_residual_j,
        });
    }
    report.bulk_amount_after_mol = inventory[species];
    report.surface_amount_after_mol = finite(surface_total.value())?;
    report.gross_adsorption_mol = finite(adsorbed.value())?;
    report.gross_desorption_mol = finite(desorbed.value())?;
    poll(&mut cancelled)?;
    bulk.replace_inventories(bulk.solvent_volume_m3(), &inventory[..count])?;
    amounts_mol.copy_from_slice(&surface);
    Ok(report)
}

struct Exchange<'a> {
    bulk: &'a AqueousSolution,
    species: usize,
    sites: &'a [Site],
    original: &'a [f64],
    bulk_amount: f64,
    duration_s: f64,
}

#[derive(Clone, Copy)]
struct Evaluation {
    bulk_amount: f64,
    residual: f64,
    adsorption: f64,
    desorption: f64,
    kinetic_residual: f64,
}

impl Exchange<'_> {
    fn evaluate(
        &self,
        amount: f64,
        next: &mut [f64],
        controls: SolubleSurfaceControls,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation> {
        poll(cancelled)?;
        if *evaluations >= controls.maximum_root_evaluations {
            return Err(SolubleSurfaceError::RootBudget {
                completed_evaluations: *evaluations,
            });
        }
        *evaluations += 1;
        let model = self
            .bulk
            .nonideal_model()
            .ok_or(SolubleSurfaceError::Invalid("explicit bulk model required"))?;
        let count = self.bulk.solute_moles().len();
        let mut concentration = [0.0; MAX_OSMOTIC_SOLUTES];
        self.bulk.concentrations_into(&mut concentration[..count])?;
        concentration[self.species] = finite(amount / self.bulk.solvent_volume_m3())?;
        let mut excess = [0.0; MAX_OSMOTIC_SOLUTES];
        model
            .log_activity_coefficients_into(&concentration[..count], &mut excess[..count])
            .map_err(OsmosisError::from)?;
        let log_activity = if amount == 0.0 {
            0.0
        } else {
            finite(
                det::ln(concentration[self.species])
                    - det::ln(model.reference_concentration_mol_m3())
                    + excess[self.species],
            )?
        };
        let mut residual = Sum::default();
        let mut adsorption = Sum::default();
        let mut desorption = Sum::default();
        let mut kinetic_residual: f64 = 0.0;
        residual.add(amount - self.bulk_amount);
        for (face, ((site, &old), output)) in
            self.sites.iter().zip(self.original).zip(next).enumerate()
        {
            poll_at(cancelled, face)?;
            let rate = if amount == 0.0 || site.desorption == 0.0 {
                0.0
            } else {
                let rate = det::exp(finite(log_activity + site.log_adsorption_prefactor)?);
                if !rate.is_finite() || rate <= 0.0 {
                    return Err(SolubleSurfaceError::Unrepresentable);
                }
                rate
            };
            let step =
                backward_euler_site(old, site.capacity, rate, site.desorption, self.duration_s)
                    .map_err(|error| SolubleSurfaceError::Site { face, error })?;
            *output = step.occupied_moles;
            residual.add(step.occupied_change_mol);
            adsorption.add(step.gross_adsorption_mol);
            desorption.add(step.gross_desorption_mol);
            kinetic_residual = kinetic_residual.max(step.kinetic_balance_residual_mol);
        }
        if kinetic_residual > controls.root_tolerance_mol {
            return Err(SolubleSurfaceError::RootAccuracy {
                residual_mol: kinetic_residual,
            });
        }
        Ok(Evaluation {
            bulk_amount: amount,
            residual: finite(residual.value())?,
            adsorption: finite(adsorption.value())?,
            desorption: finite(desorption.value())?,
            kinetic_residual,
        })
    }

    fn solve(
        &self,
        next: &mut [f64],
        controls: SolubleSurfaceControls,
        evaluations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Evaluation> {
        let initial = self.evaluate(self.bulk_amount, next, controls, evaluations, cancelled)?;
        if initial.residual.abs() <= controls.root_tolerance_mol {
            return Ok(initial);
        }
        let model = self
            .bulk
            .nonideal_model()
            .ok_or(SolubleSurfaceError::Invalid("explicit bulk model required"))?;
        let mut total = Sum::default();
        total.add(self.bulk_amount);
        for (face, &amount) in self.original.iter().enumerate() {
            poll_at(cancelled, face)?;
            total.add(amount);
        }
        let total = finite(total.value())?;
        let cap = model.maximum_concentrations_mol_m3()[self.species];
        let volume = self.bulk.solvent_volume_m3();
        let upper = if cap >= total / volume {
            total
        } else {
            let boundary = finite(cap * volume)?;
            // Keep the bracket within the source's inclusive concentration
            // support even if multiplication and division round oppositely.
            if boundary / volume > cap {
                boundary.next_down()
            } else {
                boundary
            }
        };
        let (mut low, mut high, boundary) = if initial.residual > 0.0 {
            (0.0, self.bulk_amount, 0.0)
        } else {
            (self.bulk_amount, upper, upper)
        };
        let bound = self.evaluate(boundary, next, controls, evaluations, cancelled)?;
        if bound.residual.abs() <= controls.root_tolerance_mol {
            return Ok(bound);
        }
        if bound.residual.signum() == initial.residual.signum() {
            return Err(SolubleSurfaceError::BulkSupport {
                residual_mol: bound.residual,
            });
        }
        let mut last_residual = bound.residual;
        loop {
            let middle = low + 0.5 * (high - low);
            if middle == low || middle == high {
                return Err(SolubleSurfaceError::RootAccuracy {
                    residual_mol: last_residual,
                });
            }
            let trial = self.evaluate(middle, next, controls, evaluations, cancelled)?;
            last_residual = trial.residual;
            if trial.residual.abs() <= controls.root_tolerance_mol {
                return Ok(trial);
            }
            if trial.residual > 0.0 {
                high = middle;
            } else {
                low = middle;
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
// Shared original endpoint account: a moving consumer supplies both actual
// areas without freezing the final amount on an inadmissible old geometry.
pub(in crate::capillary::surface) fn energy_change(
    bulk: &AqueousSolution,
    species: usize,
    before_bulk: &[f64],
    after_bulk: &[f64],
    laws: &[SurfactantLaw],
    before_areas: &[f64],
    after_areas: &[f64],
    before: &[f64],
    after: &[f64],
    offset: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64> {
    let model = bulk
        .nonideal_model()
        .ok_or(SolubleSurfaceError::Invalid("explicit bulk model required"))?;
    let volume = bulk.solvent_volume_m3();
    let mut change = Sum::default();
    change.add(
        model
            .reference_cancelled_inventory_free_energy_change_j(
                volume,
                before_bulk,
                volume,
                after_bulk,
            )
            .map_err(OsmosisError::from)?,
    );
    // Bulk molecules are not conserved separately: restore the original
    // reference term omitted by the source's explicitly reference-cancelled API.
    change.add(finite(
        -R * bulk.temperature_k()
            * det::ln(model.reference_concentration_mol_m3())
            * (after_bulk[species] - before_bulk[species]),
    )?);
    for (face, (((&law, &before_area), &after_area), (&old, &new))) in laws
        .iter()
        .zip(before_areas)
        .zip(after_areas)
        .zip(before.iter().zip(after))
        .enumerate()
    {
        poll_at(cancelled, face)?;
        change.add(
            SurfactantFace {
                law,
                amount_mol: old,
            }
            .energy_change(before_area, after_area, new)
            .map_err(|error| SolubleSurfaceError::Surface { face, error })?,
        );
        change.add(finite(offset * (new - old))?);
    }
    finite(change.value())
}

fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(SolubleSurfaceError::Unrepresentable)
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SolubleSurfaceError::Cancelled)
    } else {
        Ok(())
    }
}
fn poll_at(cancelled: &mut impl FnMut() -> bool, index: usize) -> Result<()> {
    if index % 64 == 0 {
        poll(cancelled)?;
    }
    Ok(())
}
fn copied(values: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<Vec<f64>> {
    let mut copy = Vec::new();
    copy.try_reserve_exact(values.len())
        .map_err(|_| SolubleSurfaceError::Allocation)?;
    for (index, &value) in values.iter().enumerate() {
        poll_at(cancelled, index)?;
        copy.push(value);
    }
    Ok(copy)
}
fn admit_energy(change_j: f64, controls: SolubleSurfaceControls) -> Result<()> {
    if change_j > controls.energy_tolerance_j {
        Err(SolubleSurfaceError::FreeEnergy { change_j })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
