//! Conservative density/inventory continuation between prepared native meshes.
//!
//! The original mesh remapper moves physical volume, mass and every `mass*q`
//! through the SAME caller-supplied overlaps. Receiving native cell volumes
//! determine density; received mass determines each specific value. The target
//! plan already owns its separately prescribed, admitted volume advector.
//! This transfers P0 material inventories, not the BDM velocity/pressure field,
//! angular momentum, kinetic energy, constitutive history or overlap geometry.

use fs_exec::Cx;
use fs_mesh::{
    ConservativeRemapError, ConservativeRemapOutcome, MAX_REMAP_SOURCE_COVERAGE_TOLERANCE,
    RemapContribution, conservative_cell_remap,
};

use super::{
    DensityTransportCell, DensityTransportError, PreparedDensityTransport, Sum, admit_cell,
    admit_storage, finite, positive_value,
};

/// Admission of one complete P0 material-state transfer.
#[derive(Clone, Copy, Debug)]
pub struct DensityRemapControls<const N: usize> {
    /// Absolute source-row fraction defect, in [0, 1e-6], as in `fs-mesh`.
    pub source_coverage_tolerance: f64,
    /// Per-target |mapped volume - native volume| / max(both volumes), in
    /// [0, 1e-6]. A target merely receiving some overlap is not full coverage.
    pub relative_target_volume_tolerance: f64,
    /// Absolute total physical-volume defect [m³], for original overlap
    /// transfer AND the actual source/target native volumes.
    pub volume_tolerance_m3: f64,
    /// Absolute mass defect [kg], including reconstructed density roundoff.
    pub mass_tolerance_kg: f64,
    /// Separate absolute defect allowance for each extensive `mass*specific`,
    /// in its original declared units; signed inventories are supported.
    pub carried_tolerance: [f64; N],
    /// Adapter-owned vector payload: reusable source scalar, retained target
    /// masses, target states, and one temporary remapped scalar vector; also
    /// the separate volume phase's remapped values and compensated coverages.
    /// The original remapper additionally enforces its own per-call cell,
    /// overlap and auxiliary-memory caps. This is not a whole-process RSS cap.
    pub maximum_storage_bytes: usize,
}

/// Physical field whose transfer or conservation allowance refused the request.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DensityRemapQuantity {
    /// Physical cell volume [m³].
    Volume,
    /// Cell mass [kg].
    Mass,
    /// `mass*specific` for this component, in its caller-declared units.
    Carried(usize),
}

/// Refusal before any target state is published or either plan is changed.
#[derive(Clone, Debug, PartialEq)]
pub enum DensityRemapError {
    /// Native density, inventory, storage or representability admission.
    Transport(DensityTransportError),
    /// Original scalar-remap refusal with the physical field retained.
    Remap {
        quantity: DensityRemapQuantity,
        error: ConservativeRemapError,
    },
    /// A target's supplied overlaps do not represent its native physical volume.
    TargetVolumeCoverage {
        target: usize,
        mapped_m3: f64,
        native_m3: f64,
        defect_m3: f64,
        relative_defect: f64,
    },
    /// Actual reconstructed-state total differs from the source beyond allowance.
    Conservation {
        quantity: DensityRemapQuantity,
        residual: f64,
        maximum: f64,
    },
    /// Cancellation during adapter work; no candidate was published.
    Cancelled,
}

impl core::fmt::Display for DensityRemapError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native density remap: {self:?}")
    }
}
impl std::error::Error for DensityRemapError {}
impl From<DensityTransportError> for DensityRemapError {
    fn from(error: DensityTransportError) -> Self {
        Self::Transport(error)
    }
}
type Result<T> = core::result::Result<T, DensityRemapError>;

/// Source inventories and ACTUAL reconstructed native target inventories.
/// All defects are measured f64 values, not repairs or geometric certificates.
/// Residuals reduce individual source/target amounts together with compensation;
/// carried residuals also retain FMA product remainders, as native transport does.
/// They can be nonzero even when separately rounded endpoint totals are equal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityRemapReport<const N: usize> {
    pub source_volume_m3: f64,
    pub target_volume_m3: f64,
    pub volume_balance_residual_m3: f64,
    pub maximum_source_coverage_defect: f64,
    pub maximum_target_volume_defect_m3: f64,
    pub maximum_relative_target_volume_defect: f64,
    pub mass_before_kg: f64,
    pub mass_after_kg: f64,
    pub mass_balance_residual_kg: f64,
    pub carried_before: [f64; N],
    pub carried_after: [f64; N],
    pub carried_balance_residual: [f64; N],
}

/// Accepted material cells ready for the target plan's original `advance_into`.
#[derive(Clone, Debug, PartialEq)]
pub struct DensityRemapOutcome<const N: usize> {
    pub cells: Vec<DensityTransportCell<N>>,
    pub report: DensityRemapReport<N>,
}

impl PreparedDensityTransport<'_> {
    /// Restart P0 density and specific inventories on THIS prepared target plan.
    ///
    /// `source` supplies actual old volumes; this plan supplies actual new
    /// volumes and its separately admitted advector. Canonical contributions
    /// are fractions of SOURCE physical volume. They must partition every
    /// source and cover every target's native volume within the stated bounds.
    /// Geometry construction/certification remains with the caller.
    ///
    /// Each scalar uses `fs_mesh::conservative_cell_remap` unchanged. Extensive
    /// inventories are formed BEFORE division by receiving mass, preserving
    /// representable carried amounts even for very small, high-specific donors.
    /// The final audit recomputes `rho*V` and `(rho*V)*q` from the returned cells.
    /// No density clipping, species renormalization or energy repair occurs.
    ///
    /// Input cells and both plans remain unchanged on every error/cancellation.
    /// A successful candidate can immediately enter this plan's density Euler
    /// or subcycled transport. No coupled Navier--Stokes checkpoint, velocity
    /// projection, angular/kinetic conservation or physical clock is inferred.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn remap_from<const N: usize>(
        &self,
        source: &PreparedDensityTransport<'_>,
        cells: &[DensityTransportCell<N>],
        contributions: &[RemapContribution],
        controls: DensityRemapControls<N>,
        cx: &Cx<'_>,
    ) -> Result<DensityRemapOutcome<N>> {
        checkpoint(cx)?;
        let source_volumes = source.cell_volumes_m3();
        let target_volumes = self.cell_volumes_m3();
        if cells.len() != source_volumes.len() {
            return Err(DensityTransportError::Invalid("remap source dimensions differ").into());
        }
        admit_controls(&controls)?;
        admit_storage(
            &[
                (target_volumes.len(), size_of::<f64>()),
                (target_volumes.len(), size_of::<Sum>()),
            ],
            controls.maximum_storage_bytes,
        )?;
        // At most one scalar result coexists with source scratch, target mass
        // and target states. The volume result is dropped before those allocate.
        admit_storage(
            &[
                (cells.len(), size_of::<f64>()),
                (target_volumes.len(), size_of::<f64>()),
                (target_volumes.len(), size_of::<f64>()),
                (target_volumes.len(), size_of::<DensityTransportCell<N>>()),
            ],
            controls.maximum_storage_bytes,
        )?;
        let transfer = |values: &[f64], quantity, tolerance| {
            conservative_cell_remap(
                values,
                target_volumes.len(),
                contributions,
                controls.source_coverage_tolerance,
                tolerance,
                cx,
            )
            .map_err(|error| DensityRemapError::Remap { quantity, error })
        };

        let mut report = {
            let volumes = transfer(
                source_volumes,
                DensityRemapQuantity::Volume,
                controls.volume_tolerance_m3,
            )?;
            let mut total = Sum::default();
            let mut volume_balance = Sum::default();
            for &volume in source_volumes {
                checkpoint(cx)?;
                volume_balance.add(-volume)?;
            }
            let mut coverage = reserve(target_volumes.len())?;
            for &native in target_volumes {
                checkpoint(cx)?;
                let mut difference = Sum::default();
                difference.add(-native)?;
                coverage.push(difference);
            }
            // The original owner has validated every index/fraction. Preserve
            // small local coverage defects that rounding its target value can
            // hide; these are audits, not a second set of remapped values.
            for contribution in contributions {
                checkpoint(cx)?;
                coverage[contribution.target()].add(finite(
                    source_volumes[contribution.source()] * contribution.fraction(),
                    "remap overlap volume",
                )?)?;
            }
            let mut maximum_defect: f64 = 0.0;
            let mut maximum_relative_defect: f64 = 0.0;
            for (target, (&mapped, &native)) in volumes
                .target_values()
                .iter()
                .zip(target_volumes)
                .enumerate()
            {
                checkpoint(cx)?;
                positive_value(mapped, "remap target overlap volume")?;
                let defect = coverage[target].value()?.abs();
                let scale = mapped.max(native);
                let relative = defect / scale;
                // Multiplying the allowance also retains a nonzero defect
                // when a relative quotient would underflow at zero tolerance.
                if defect > controls.relative_target_volume_tolerance * scale {
                    return Err(DensityRemapError::TargetVolumeCoverage {
                        target,
                        mapped_m3: mapped,
                        native_m3: native,
                        defect_m3: defect,
                        relative_defect: relative,
                    });
                }
                maximum_defect = maximum_defect.max(defect);
                maximum_relative_defect = maximum_relative_defect.max(relative);
                total.add(native)?;
                volume_balance.add(native)?;
            }
            let before = volumes.report().source_total();
            let after = total.value()?;
            DensityRemapReport {
                source_volume_m3: before,
                target_volume_m3: after,
                volume_balance_residual_m3: balance(
                    DensityRemapQuantity::Volume,
                    volume_balance,
                    controls.volume_tolerance_m3,
                )?,
                maximum_source_coverage_defect: volumes.report().maximum_source_coverage_defect(),
                maximum_target_volume_defect_m3: maximum_defect,
                maximum_relative_target_volume_defect: maximum_relative_defect,
                mass_before_kg: 0.0,
                mass_after_kg: 0.0,
                mass_balance_residual_kg: 0.0,
                carried_before: [0.0; N],
                carried_after: [0.0; N],
                carried_balance_residual: [0.0; N],
            }
        };

        let mut values = reserve(cells.len())?;
        let mut mass_balance = Sum::default();
        for (cell, &volume) in cells.iter().zip(source_volumes) {
            checkpoint(cx)?;
            admit_cell(cell)?;
            let mass = positive_value(cell.density_kg_m3 * volume, "remap source mass")?;
            mass_balance.add(-mass)?;
            values.push(mass);
        }
        let mass = transfer(
            &values,
            DensityRemapQuantity::Mass,
            controls.mass_tolerance_kg,
        )?;
        report.mass_before_kg = mass.report().source_total();
        let masses = mass.into_target_values();
        let mut output = reserve(target_volumes.len())?;
        let mut actual_mass = Sum::default();
        for (&mass, &volume) in masses.iter().zip(target_volumes) {
            checkpoint(cx)?;
            positive_value(mass, "remap received mass")?;
            let density = positive_value(mass / volume, "remap reconstructed density")?;
            let represented = positive_value(density * volume, "remap represented mass")?;
            actual_mass.add(represented)?;
            mass_balance.add(represented)?;
            output.push(DensityTransportCell {
                density_kg_m3: density,
                specific: [0.0; N],
            });
        }
        report.mass_after_kg = actual_mass.value()?;
        report.mass_balance_residual_kg = balance(
            DensityRemapQuantity::Mass,
            mass_balance,
            controls.mass_tolerance_kg,
        )?;

        for component in 0..N {
            checkpoint(cx)?;
            let mut carried_balance = Sum::default();
            for (index, value) in values.iter_mut().enumerate() {
                checkpoint(cx)?;
                let mass = positive_value(
                    cells[index].density_kg_m3 * source_volumes[index],
                    "remap source mass",
                )?;
                *value = finite(
                    mass * cells[index].specific[component],
                    "remap source inventory",
                )?;
                carried_balance.add(-*value)?;
                // Audit the product's low part too, as native advance_into
                // does: a large common specific datum must not hide a small
                // reconstruction defect behind equal rounded inventories.
                carried_balance.add(-mass.mul_add(cells[index].specific[component], -*value))?;
            }
            let quantity = DensityRemapQuantity::Carried(component);
            let carried: ConservativeRemapOutcome =
                transfer(&values, quantity, controls.carried_tolerance[component])?;
            report.carried_before[component] = carried.report().source_total();
            let mut actual = Sum::default();
            for (index, cell) in output.iter_mut().enumerate() {
                checkpoint(cx)?;
                cell.specific[component] = finite(
                    carried.target_values()[index] / masses[index],
                    "remap reconstructed specific inventory",
                )?;
                let mass = positive_value(
                    cell.density_kg_m3 * target_volumes[index],
                    "remap represented mass",
                )?;
                let represented = finite(
                    mass * cell.specific[component],
                    "remap represented inventory",
                )?;
                actual.add(represented)?;
                carried_balance.add(represented)?;
                carried_balance.add(mass.mul_add(cell.specific[component], -represented))?;
            }
            report.carried_after[component] = actual.value()?;
            report.carried_balance_residual[component] = balance(
                quantity,
                carried_balance,
                controls.carried_tolerance[component],
            )?;
        }
        checkpoint(cx)?;
        Ok(DensityRemapOutcome {
            cells: output,
            report,
        })
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint().map_err(|_| DensityRemapError::Cancelled)
}

fn reserve<T>(count: usize) -> Result<Vec<T>> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| DensityTransportError::Allocation)?;
    Ok(values)
}

fn admit_controls<const N: usize>(controls: &DensityRemapControls<N>) -> Result<()> {
    if [
        controls.source_coverage_tolerance,
        controls.relative_target_volume_tolerance,
    ]
    .iter()
    .any(|value| !value.is_finite() || *value < 0.0 || *value > MAX_REMAP_SOURCE_COVERAGE_TOLERANCE)
        || [controls.volume_tolerance_m3, controls.mass_tolerance_kg]
            .iter()
            .chain(&controls.carried_tolerance)
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(DensityTransportError::Invalid(
            "finite nonnegative remap allowances required; coverage at most 1e-6",
        )
        .into());
    }
    Ok(())
}

fn balance(quantity: DensityRemapQuantity, difference: Sum, maximum: f64) -> Result<f64> {
    let residual = difference.value()?;
    if residual.abs() > maximum {
        return Err(DensityRemapError::Conservation {
            quantity,
            residual,
            maximum,
        });
    }
    Ok(residual)
}
