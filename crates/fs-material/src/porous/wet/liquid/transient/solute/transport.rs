//! Conservative passive neutral-solute transport on one accepted water interval.
//! Amounts and flux integrals are per unit bulk cross-sectional area.

use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};

pub(super) mod nonlinear;

/// Independent solute boundary; the hydraulic owner has already selected water flux.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiluteSoluteBoundary<const N: usize> {
    /// No molecular exchange; zero or outward water flow is admitted.
    Outflow,
    /// No molecular exchange. This concentration supplies water entering the
    /// column; outward water carries the actual interior concentration.
    InflowConcentration([f64; N]),
    /// Fixed external concentration with independent film resistances [s/m].
    /// j_out=(c_surface-c_external)/R. Zero R prescribes the actual face value;
    /// cell D=0 still blocks diffusion. Incoming water uses this same bath.
    Reservoir {
        concentration_mol_m3: [f64; N],
        film_resistance_s_m: [f64; N],
    },
}

/// Numerical residual and conservation budgets; none is a temporal error bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiluteSoluteControls<const N: usize> {
    /// Per-species allowance on EVERY original cell amount equation [mol/m²].
    /// The sorbing entry also applies it to initial equilibrium consistency and
    /// represented bound exchange; its transport equation includes both pools.
    pub absolute_equation_tolerance_mol_m2: [f64; N],
    /// Additional allowance on actual storage change or incident face transfer.
    /// Background inventory does not enlarge this scale. In [0,1e-3].
    pub relative_equation_tolerance: f64,
    /// Fixed whole-column constituent balance allowances [mol/m²].
    pub amount_balance_tolerance_mol_m2: [f64; N],
    /// Per-cell AND whole-column accepted hydraulic residual allowance [m].
    pub water_balance_tolerance_m: f64,
    /// Bound on cells * species, independently of physical timestep.
    pub maximum_cell_species: usize,
    /// Conservative additional numerical payload, including output arrays and
    /// original/final water views; excludes existing carriers/inputs and allocator
    /// metadata. The native tridiagonal factors/scratch are fully included.
    pub maximum_storage_bytes: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DiluteSoluteReport<const N: usize> {
    pub interval_s: [f64; 2],
    pub linear_solves: usize,
    pub cell_species: usize,
    pub workspace_bytes: usize,
    pub maximum_water_balance_defect_m: f64,
    pub amount_before_mol_m2: [f64; N],
    pub amount_after_mol_m2: [f64; N],
    /// Compensated sum of actual individual cell amount changes.
    pub amount_change_mol_m2: [f64; N],
    pub boundary_amount_into_column_mol_m2: [f64; N],
    pub amount_balance_defect_mol_m2: [f64; N],
    pub maximum_equation_defect_mol_m2: [f64; N],
    /// Signed left-to-right integrals at all cells+1 faces [mol/m²].
    /// These are distinct mechanisms; neither should be applied a second time.
    pub advective_face_amount_mol_m2: Vec<[f64; N]>,
    pub diffusive_face_amount_mol_m2: Vec<[f64; N]>,
}

/// Private numerical proposal; the borrowed hydraulic candidate is unchanged.
#[derive(Clone, Debug, PartialEq)]
pub struct DiluteSoluteStep<const N: usize> {
    pub amount_mol_m2: Vec<[f64; N]>,
    pub concentration_mol_m3: Vec<[f64; N]>,
    pub report: DiluteSoluteReport<N>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DiluteSoluteError {
    Invalid(&'static str),
    Workspace {
        required_bytes: usize,
        maximum_bytes: usize,
    },
    WorkBudget {
        required_cell_species: usize,
        maximum_cell_species: usize,
    },
    Linear(TridiagonalError),
    WaterBalance {
        cell: Option<usize>,
        residual_m: f64,
    },
    Equation {
        cell: usize,
        species: usize,
        residual_mol_m2: f64,
        allowance_mol_m2: f64,
    },
    Sorption {
        cell: usize,
        species: usize,
        residual_mol_m2: f64,
        allowance_mol_m2: f64,
    },
    Nonlinear {
        species: usize,
        linear_solves: usize,
        maximum_residual_mol_m2: f64,
    },
    NonlinearWork {
        required_linear_solves: usize,
        maximum_linear_solves: usize,
    },
    Balance {
        species: usize,
        residual_mol_m2: f64,
        allowance_mol_m2: f64,
    },
    Cancelled,
}
impl core::fmt::Display for DiluteSoluteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "porous dilute-solute transport: {self:?}")
    }
}
impl std::error::Error for DiluteSoluteError {}
impl From<TridiagonalError> for DiluteSoluteError {
    fn from(value: TridiagonalError) -> Self {
        if value == TridiagonalError::Cancelled {
            Self::Cancelled
        } else {
            Self::Linear(value)
        }
    }
}
pub(super) type Result<T> = core::result::Result<T, DiluteSoluteError>;
pub(super) fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
pub(super) fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
pub(super) fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(DiluteSoluteError::Invalid("nonfinite solute arithmetic"))
    }
}
pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(DiluteSoluteError::Cancelled)
    } else {
        Ok(())
    }
}
pub(super) fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64> {
    let (mut total, mut correction) = (0.0, 0.0);
    for x in values {
        let y = x - correction;
        let next = total + y;
        correction = (next - total) - y;
        total = next;
    }
    finite(total)
}

pub(super) fn admit_extent<const N: usize>(
    n: usize,
    controls: DiluteSoluteControls<N>,
) -> Result<(usize, usize)> {
    if n == 0
        || n > 4096
        || N == 0
        || N > 32
        || !controls
            .absolute_equation_tolerance_mol_m2
            .iter()
            .all(|&v| positive(v))
        || !controls
            .amount_balance_tolerance_mol_m2
            .iter()
            .all(|&v| positive(v))
        || !positive(controls.water_balance_tolerance_m)
        || !nonnegative(controls.relative_equation_tolerance)
        || controls.relative_equation_tolerance > 1e-3
    {
        return Err(DiluteSoluteError::Invalid(
            "bounded wet column, positive budgets and relative tolerance in [0,1e-3] required",
        ));
    }
    let work = n * N;
    if work > controls.maximum_cell_species {
        return Err(DiluteSoluteError::WorkBudget {
            required_cell_species: work,
            maximum_cell_species: controls.maximum_cell_species,
        });
    }
    // Four field/face arrays, adapter water/geometry views, one reusable
    // matrix/RHS and all five native factor/scratch arrays, plus bounded scalar
    // and report payload. Vector capacities are requested exactly.
    let bytes = (4 * (n + 1) * N + 15 * (n + 1) + 16 * N) * core::mem::size_of::<f64>()
        + n * core::mem::size_of::<bool>()
        + 1024;
    if bytes > controls.maximum_storage_bytes {
        return Err(DiluteSoluteError::Workspace {
            required_bytes: bytes,
            maximum_bytes: controls.maximum_storage_bytes,
        });
    }
    Ok((work, bytes))
}

impl<const N: usize> DiluteSoluteBoundary<N> {
    fn values(self) -> Option<[f64; N]> {
        match self {
            Self::Outflow => None,
            Self::InflowConcentration(c) => Some(c),
            Self::Reservoir {
                concentration_mol_m3: c,
                ..
            } => Some(c),
        }
    }
    fn admit(self, inflow: bool) -> Result<()> {
        if self
            .values()
            .is_some_and(|c| !c.iter().all(|&v| nonnegative(v)))
        {
            return Err(DiluteSoluteError::Invalid(
                "finite nonnegative boundary concentrations required",
            ));
        }
        if let Self::Reservoir {
            film_resistance_s_m,
            ..
        } = self
            && !film_resistance_s_m.iter().all(|&v| nonnegative(v))
        {
            return Err(DiluteSoluteError::Invalid(
                "finite nonnegative diffusion-film resistances required",
            ));
        }
        if inflow && self.values().is_none() {
            return Err(DiluteSoluteError::Invalid(
                "water entering the column requires explicit solute concentrations",
            ));
        }
        Ok(())
    }
}

/// Reciprocal of two positive/nonnegative resistances, avoiding sum overflow.
fn series(left: f64, right: f64) -> Result<f64> {
    if !positive(left) || !nonnegative(right) {
        return Err(DiluteSoluteError::Invalid(
            "wet half-cell diffusion resistance is unrepresentable",
        ));
    }
    let largest = left.max(right);
    let g = (1.0 / (1.0 + left.min(right) / largest)) / largest;
    if positive(g) {
        Ok(g)
    } else {
        Err(DiluteSoluteError::Invalid(
            "positive face conductance is unrepresentable",
        ))
    }
}

pub(super) fn integrated_flux(coefficient: f64, drive: f64) -> Result<f64> {
    let value = finite(coefficient * drive)?;
    if coefficient != 0.0 && drive != 0.0 && value == 0.0 {
        return Err(DiluteSoluteError::Invalid(
            "nonzero solute transfer is unrepresentable",
        ));
    }
    Ok(value)
}

/// Optional local immobile storage. It never changes the water geometry or
/// any face coefficient. The zero-sized free-solute policy allocates nothing.
pub(super) trait ImmobileStorage<const N: usize> {
    fn additional_bytes(&self) -> usize {
        0
    }
    fn admit(
        &mut self,
        _: &[[f64; N]],
        _: &[f64],
        _: &[f64],
        _: DiluteSoluteControls<N>,
        _: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        Ok(())
    }
    /// Add only to the concentration-system storage diagonal and amount RHS.
    fn terms(&self, _: usize, _: usize) -> Result<(f64, f64)> {
        Ok((0.0, 0.0))
    }
    /// Solve the assembled ORIGINAL aqueous stencil. Nonlinear immobile
    /// storage can change its local diagonal, never its physical face terms.
    #[allow(clippy::too_many_arguments)]
    fn solve_species(
        &mut self,
        _: usize,
        lower: &[f64],
        diagonal: &[f64],
        upper: &[f64],
        rhs: &mut [f64],
        solver: &mut TridiagonalWorkspace,
        _: DiluteSoluteControls<N>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize> {
        solver.factor(lower, diagonal, upper, &mut *cancelled)?;
        solver.solve_factored_in_place(rhs, &mut *cancelled)?;
        Ok(1)
    }
    fn finish_cell(
        &mut self,
        _: usize,
        _: usize,
        _: f64,
        _: DiluteSoluteControls<N>,
    ) -> Result<()> {
        Ok(())
    }
    /// Original amount, represented final amount, represented amount change.
    fn amounts(&self, _: usize, _: usize) -> [f64; 3] {
        [0.0; 3]
    }
}
impl<const N: usize> ImmobileStorage<N> for () {}

#[cfg(test)]
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn transport<const N: usize>(
    interval_s: [f64; 2],
    duration_s: f64,
    before_water_m: &[f64],
    after_water_m: &[f64],
    widths_m: &[f64],
    water_face_m: &[f64],
    initial: &[[f64; N]],
    diffusivity_m2_s: &[[f64; N]],
    tortuosity: &[f64],
    boundaries: [DiluteSoluteBoundary<N>; 2],
    controls: DiluteSoluteControls<N>,
    cancelled: impl FnMut() -> bool,
) -> Result<DiluteSoluteStep<N>> {
    transport_with_storage(
        interval_s,
        duration_s,
        before_water_m,
        after_water_m,
        widths_m,
        water_face_m,
        initial,
        diffusivity_m2_s,
        tortuosity,
        boundaries,
        controls,
        &mut (),
        cancelled,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn transport_with_storage<const N: usize>(
    interval_s: [f64; 2],
    duration_s: f64,
    before_water_m: &[f64],
    after_water_m: &[f64],
    widths_m: &[f64],
    water_face_m: &[f64],
    initial: &[[f64; N]],
    diffusivity_m2_s: &[[f64; N]],
    tortuosity: &[f64],
    boundaries: [DiluteSoluteBoundary<N>; 2],
    controls: DiluteSoluteControls<N>,
    storage: &mut impl ImmobileStorage<N>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DiluteSoluteStep<N>> {
    poll(&mut cancelled)?;
    let n = initial.len();
    let (cell_species, original_bytes) = admit_extent(n, controls)?;
    let workspace_bytes = original_bytes
        .checked_add(storage.additional_bytes())
        .ok_or(DiluteSoluteError::Invalid(
            "sorption storage extent overflow",
        ))?;
    if workspace_bytes > controls.maximum_storage_bytes {
        return Err(DiluteSoluteError::Workspace {
            required_bytes: workspace_bytes,
            maximum_bytes: controls.maximum_storage_bytes,
        });
    }
    if before_water_m.len() != n
        || after_water_m.len() != n
        || widths_m.len() != n
        || water_face_m.len() != n + 1
        || diffusivity_m2_s.len() != n
        || tortuosity.len() != n
        || !positive(duration_s)
        || !interval_s.iter().all(|v| v.is_finite())
        || interval_s[1] <= interval_s[0]
    {
        return Err(DiluteSoluteError::Invalid(
            "matching transport extents and advancing interval required",
        ));
    }
    boundaries[0].admit(water_face_m[0] > 0.0)?;
    boundaries[1].admit(water_face_m[n] < 0.0)?;
    if !water_face_m.iter().all(|v| v.is_finite()) {
        return Err(DiluteSoluteError::Invalid(
            "finite original hydraulic transfers required",
        ));
    }
    let mut maximum_water_balance_defect_m = 0.0_f64;
    for i in 0..n {
        if i % 128 == 0 {
            poll(&mut cancelled)?;
        }
        if !positive(before_water_m[i])
            || !positive(after_water_m[i])
            || !positive(widths_m[i])
            || after_water_m[i] > widths_m[i]
            || before_water_m[i] > widths_m[i]
            || !tortuosity[i].is_finite()
            || tortuosity[i] < 1.0
            || !initial[i].iter().all(|&v| nonnegative(v))
            || !diffusivity_m2_s[i].iter().all(|&v| nonnegative(v))
        {
            return Err(DiluteSoluteError::Invalid(
                "strictly wet storage, valid porosity, nonnegative amounts/D and transport tortuosity >=1 required",
            ));
        }
        let residual = finite(
            (after_water_m[i] - before_water_m[i]) - (water_face_m[i] - water_face_m[i + 1]),
        )?;
        maximum_water_balance_defect_m = maximum_water_balance_defect_m.max(residual.abs());
        if residual.abs() > controls.water_balance_tolerance_m {
            return Err(DiluteSoluteError::WaterBalance {
                cell: Some(i),
                residual_m: residual,
            });
        }
    }
    let water_residual = finite(
        sum((0..n).map(|i| after_water_m[i] - before_water_m[i]))?
            - (water_face_m[0] - water_face_m[n]),
    )?;
    maximum_water_balance_defect_m = maximum_water_balance_defect_m.max(water_residual.abs());
    if water_residual.abs() > controls.water_balance_tolerance_m {
        return Err(DiluteSoluteError::WaterBalance {
            cell: None,
            residual_m: water_residual,
        });
    }
    storage.admit(
        initial,
        before_water_m,
        after_water_m,
        controls,
        &mut cancelled,
    )?;
    let mut step = DiluteSoluteStep {
        amount_mol_m2: vec![[0.0; N]; n],
        concentration_mol_m3: vec![[0.0; N]; n],
        report: DiluteSoluteReport {
            interval_s,
            linear_solves: 0,
            cell_species,
            workspace_bytes,
            maximum_water_balance_defect_m,
            amount_before_mol_m2: [0.0; N],
            amount_after_mol_m2: [0.0; N],
            amount_change_mol_m2: [0.0; N],
            boundary_amount_into_column_mol_m2: [0.0; N],
            amount_balance_defect_mol_m2: [0.0; N],
            maximum_equation_defect_mol_m2: [0.0; N],
            advective_face_amount_mol_m2: vec![[0.0; N]; n + 1],
            diffusive_face_amount_mol_m2: vec![[0.0; N]; n + 1],
        },
    };
    let mut lower = vec![0.0; n - 1];
    let mut diagonal = vec![0.0; n];
    let mut upper = vec![0.0; n - 1];
    let mut rhs = vec![0.0; n];
    let mut resistance = vec![0.0; n];
    let mut conductance_dt = vec![0.0; n + 1];
    let mut solver = TridiagonalWorkspace::new(n)?;
    for k in 0..N {
        poll(&mut cancelled)?;
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            let d = diffusivity_m2_s[i][k];
            if d == 0.0 {
                resistance[i] = f64::INFINITY;
                continue;
            }
            let coefficient = (after_water_m[i] / widths_m[i]) * d / tortuosity[i];
            resistance[i] = (0.5 * widths_m[i]) / coefficient;
            if !positive(coefficient) || !positive(resistance[i]) {
                return Err(DiluteSoluteError::Invalid(
                    "positive wet-cell diffusion coefficient is unrepresentable",
                ));
            }
        }
        conductance_dt.fill(0.0);
        for face in 1..n {
            if face % 128 == 0 {
                poll(&mut cancelled)?;
            }
            if resistance[face - 1].is_finite() && resistance[face].is_finite() {
                conductance_dt[face] =
                    integrated_flux(duration_s, series(resistance[face - 1], resistance[face])?)?;
            }
        }
        for (side, face, cell) in [(0, 0, 0), (1, n, n - 1)] {
            if let DiluteSoluteBoundary::Reservoir {
                film_resistance_s_m,
                ..
            } = boundaries[side]
                && resistance[cell].is_finite()
            {
                conductance_dt[face] = integrated_flux(
                    duration_s,
                    series(resistance[cell], film_resistance_s_m[k])?,
                )?;
            }
        }
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            let (left, right) = (water_face_m[i], water_face_m[i + 1]);
            let (capacity, release) = storage.terms(i, k)?;
            diagonal[i] = finite(
                after_water_m[i]
                    + (-left).max(0.0)
                    + right.max(0.0)
                    + conductance_dt[i]
                    + conductance_dt[i + 1]
                    + capacity,
            )?;
            rhs[i] = finite(initial[i][k] + release)?;
            if i > 0 {
                lower[i - 1] = finite(-left.max(0.0) - conductance_dt[i])?;
            }
            if i + 1 < n {
                upper[i] = finite(-(-right).max(0.0) - conductance_dt[i + 1])?;
            }
        }
        for (side, face, cell, sign) in [(0, 0, 0, 1.0), (1, n, n - 1, -1.0)] {
            if let Some(c) = boundaries[side].values() {
                let coefficient =
                    finite((sign * water_face_m[face]).max(0.0) + conductance_dt[face])?;
                rhs[cell] = finite(rhs[cell] + integrated_flux(coefficient, c[k])?)?;
            }
        }
        step.report.linear_solves += storage.solve_species(
            k,
            &lower,
            &diagonal,
            &upper,
            &mut rhs,
            &mut solver,
            controls,
            &mut cancelled,
        )?;
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            let c = rhs[i];
            if !nonnegative(c) || (initial[i][k] > 0.0 && c == 0.0) {
                return Err(DiluteSoluteError::Invalid(
                    "positive implicit concentration is unrepresentable; no clipping is applied",
                ));
            }
            step.concentration_mol_m3[i][k] = c;
            step.amount_mol_m2[i][k] = integrated_flux(after_water_m[i], c)?;
            storage.finish_cell(i, k, c, controls)?;
        }
        for face in 0..=n {
            if face % 128 == 0 {
                poll(&mut cancelled)?;
            }
            let left = if face == 0 {
                boundaries[0].values().map_or(0.0, |c| c[k])
            } else {
                rhs[face - 1]
            };
            let right = if face == n {
                boundaries[1].values().map_or(0.0, |c| c[k])
            } else {
                rhs[face]
            };
            let flow = water_face_m[face];
            step.report.advective_face_amount_mol_m2[face][k] =
                integrated_flux(flow, if flow >= 0.0 { left } else { right })?;
            step.report.diffusive_face_amount_mol_m2[face][k] =
                integrated_flux(conductance_dt[face], finite(left - right)?)?;
        }
        let flux = |face: usize| {
            finite(
                step.report.advective_face_amount_mol_m2[face][k]
                    + step.report.diffusive_face_amount_mol_m2[face][k],
            )
        };
        let mut maximum_defect = 0.0_f64;
        for i in 0..n {
            if i % 128 == 0 {
                poll(&mut cancelled)?;
            }
            let change =
                finite(step.amount_mol_m2[i][k] - initial[i][k] + storage.amounts(i, k)[2])?;
            let (left, right) = (flux(i)?, flux(i + 1)?);
            let residual = finite(change - (left - right))?;
            maximum_defect = maximum_defect.max(residual.abs());
            let allowance = finite(
                controls.absolute_equation_tolerance_mol_m2[k]
                    + controls.relative_equation_tolerance
                        * change.abs().max(left.abs()).max(right.abs()),
            )?;
            if residual.abs() > allowance {
                return Err(DiluteSoluteError::Equation {
                    cell: i,
                    species: k,
                    residual_mol_m2: residual,
                    allowance_mol_m2: allowance,
                });
            }
        }
        let boundary = finite(flux(0)? - flux(n)?)?;
        // Preserve the original free-solute compensated reduction order.
        let change = if storage.additional_bytes() == 0 {
            sum((0..n).map(|i| step.amount_mol_m2[i][k] - initial[i][k]))?
        } else {
            sum((0..n).flat_map(|i| {
                [
                    step.amount_mol_m2[i][k] - initial[i][k],
                    storage.amounts(i, k)[2],
                ]
            }))?
        };
        let residual = finite(change - boundary)?;
        if residual.abs() > controls.amount_balance_tolerance_mol_m2[k] {
            return Err(DiluteSoluteError::Balance {
                species: k,
                residual_mol_m2: residual,
                allowance_mol_m2: controls.amount_balance_tolerance_mol_m2[k],
            });
        }
        step.report.amount_before_mol_m2[k] = if storage.additional_bytes() == 0 {
            sum(initial.iter().map(|a| a[k]))?
        } else {
            sum((0..n).flat_map(|i| [initial[i][k], storage.amounts(i, k)[0]]))?
        };
        step.report.amount_after_mol_m2[k] = if storage.additional_bytes() == 0 {
            sum(step.amount_mol_m2.iter().map(|a| a[k]))?
        } else {
            sum((0..n).flat_map(|i| [step.amount_mol_m2[i][k], storage.amounts(i, k)[1]]))?
        };
        step.report.amount_change_mol_m2[k] = change;
        step.report.boundary_amount_into_column_mol_m2[k] = boundary;
        step.report.amount_balance_defect_mol_m2[k] = residual;
        step.report.maximum_equation_defect_mol_m2[k] = maximum_defect;
    }
    poll(&mut cancelled)?;
    Ok(step)
}

#[cfg(test)]
mod tests;
