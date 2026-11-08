//! Absorbing exit media deposited on the original layered thermal timeline.
//!
//! A finite planar thermal column occupies a declared depth interval of one
//! original homogeneous exit half-space. Every cell receives its exact
//! whole-periodic-cell modal absorption integral. No sampled heat plane,
//! fictitious optical interface, extra FFT, or depth quadrature is introduced.
//!
//! The thermal cross section is the complete periodic optical cell. It retains
//! one native temperature per depth cell; lateral temperature gradients and
//! lateral heat transfer are not resolved. Source factors remain coherent
//! internally and add only as second moments. Expected heat drives a single
//! nonlinear thermal trajectory, not an ensemble-mean temperature.
//!
//! The supplied source and exit coefficients are HELD. Cycle/ensemble-mean
//! power is held over each thermal interval; carrier transients, beats,
//! constitutive refresh, electromagnetic storage change and parameter work
//! are not inferred. Native nonlinear capacity, conduction, events, authored
//! sources and energy acceptance retain their original owners. These absorbed
//! joules must not also be authored as a second thermal source.

use super::{
    BoundLayeredThermalBody, CellId, DepositionSimulation, LayerConductionConfig,
    LayeredThermalSolveControls, LayeredThermalSolveReport, LayeredThermalStep,
    LayeredThermalWorkspace, LayeredWaveHeatingError, admit_context, checkpoint, poll,
};
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumCoherence, AngularSpectrumError, AngularSpectrumGrid, PlanarAngularSpectrum,
    },
    layers::PlanarDielectricLayers,
};
use fs_math::dd::Dd;

/// Coherent ordinary/extraordinary absorption on the same thermal-cell owner.
pub mod crystal;

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::PlanarAngularSpectrum<'_, '_> {}
    impl Sealed for super::AngularSpectrumCoherence<'_, '_, '_> {}
}

/// Original monochromatic sources with a native analytic exit-loss integral.
/// Sealed implementations retain their own source, geometry and coherence;
/// arbitrary power callbacks cannot bypass material or field admission.
pub trait SubstrateAbsorptionSource: sealed::Sealed {
    #[doc(hidden)]
    fn original_stack(&self) -> &PlanarDielectricLayers<'_>;
    #[doc(hidden)]
    fn original_grid(&self) -> AngularSpectrumGrid;
    #[doc(hidden)]
    fn mode_factor_count(&self) -> Result<usize, LayeredWaveHeatingError>;
    #[doc(hidden)]
    fn interface_power_w(&self) -> f64;
    #[doc(hidden)]
    fn absorbed_between_w(
        &self,
        start_m: f64,
        end_m: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError>;
}

impl SubstrateAbsorptionSource for PlanarAngularSpectrum<'_, '_> {
    fn original_stack(&self) -> &PlanarDielectricLayers<'_> {
        self.source()
    }
    fn original_grid(&self) -> AngularSpectrumGrid {
        self.grid()
    }
    fn mode_factor_count(&self) -> Result<usize, LayeredWaveHeatingError> {
        grid_points(self.grid())
    }
    fn interface_power_w(&self) -> f64 {
        self.power().transmitted_w
    }
    fn absorbed_between_w(
        &self,
        start_m: f64,
        end_m: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        self.transmitted_absorbed_power_w(start_m, end_m, cancelled)
    }
}

impl SubstrateAbsorptionSource for AngularSpectrumCoherence<'_, '_, '_> {
    fn original_stack(&self) -> &PlanarDielectricLayers<'_> {
        self.source()
    }
    fn original_grid(&self) -> AngularSpectrumGrid {
        self.grid()
    }
    fn mode_factor_count(&self) -> Result<usize, LayeredWaveHeatingError> {
        grid_points(self.grid())?
            .checked_mul(self.factor_count())
            .ok_or(invalid("substrate factor work overflows"))
    }
    fn interface_power_w(&self) -> f64 {
        self.power().transmitted_w
    }
    fn absorbed_between_w(
        &self,
        start_m: f64,
        end_m: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        self.transmitted_absorbed_power_w(start_m, end_m, cancelled)
    }
}

/// Preparation and one result's logical adapter storage. Borrowed optical and
/// thermal owners, native thermal state/workspace and allocator overhead are
/// separate. Each source owner documents its mode/cell work, including boundary
/// flux queries.
#[derive(Clone, Copy, Debug)]
pub struct SubstrateWaveHeatingLimits {
    pub maximum_cells: usize,
    pub maximum_mode_cell_pairs: usize,
    pub maximum_workspace_bytes: usize,
}

/// One original physical control volume and its positive, integrated heat.
/// Depth follows the source owner: from the final optical face for an exit
/// substrate, or from the entrance face for a finite crystal film.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SubstrateCellPower {
    pub cell_id: CellId,
    pub depth_m: [f64; 2],
    pub absorbed_w: f64,
}

/// Normal power through the two boundaries of the REPRESENTED thermal body.
/// Absorption outside that body's depth interval is excluded. The substrate
/// owner excludes the preceding films; a crystal-film owner includes only its
/// represented film cells. Residual includes floating-point evaluation;
/// it is reported, never redistributed into cell heat.
/// Entering and leaving power are nonnegative sums of inward and outward net
/// flow across both faces. Each face uses the total field's normal Poynting
/// flux, including coherent interference between its incident and scattered waves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SubstrateWavePower {
    pub entering_w: f64,
    pub leaving_w: f64,
    pub absorbed_w: f64,
    pub residual_w: f64,
}

/// Joules over the actual native thermal interval. Only absorbed_j is supplied
/// to the thermal owner, and it is the sum of the actual per-cell additions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SubstrateWaveEnergy {
    pub entering_j: f64,
    pub leaving_j: f64,
    pub absorbed_j: f64,
    pub residual_j: f64,
}

/// Cached, source-bound analytic watts in original thermal-cell order.
/// Preparation costs O(cells * modes * factors); repeated thermal intervals
/// reuse these watts with O(cells) work and no electromagnetic reevaluation.
pub struct PreparedSubstrateWaveHeating<'p, 'b, 'a, S: SubstrateAbsorptionSource + ?Sized> {
    source: &'p S,
    thermal: &'b BoundLayeredThermalBody<'a>,
    cells: Vec<SubstrateCellPower>,
    power: SubstrateWavePower,
    mode_cell_pairs: usize,
    workspace_bytes: usize,
}

/// Native thermal candidate plus the distinct optical heat account. All
/// refusal/cancellation paths preserve the supplied initial simulation.
pub struct SubstrateWaveStep<'t> {
    pub thermal: DepositionSimulation<'t>,
    pub heat_report: LayeredThermalSolveReport,
    pub optical: SubstrateWaveEnergy,
    pub cell_energy_j: Vec<(CellId, f64)>,
}

fn invalid(message: &'static str) -> LayeredWaveHeatingError {
    LayeredWaveHeatingError::Invalid(message)
}

fn grid_points(grid: AngularSpectrumGrid) -> Result<usize, LayeredWaveHeatingError> {
    grid.width
        .checked_mul(grid.height)
        .ok_or(invalid("substrate grid work overflows"))
}

fn reserve<T>(n: usize) -> Result<Vec<T>, LayeredWaveHeatingError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| invalid("substrate adapter allocation failed"))?;
    Ok(values)
}

fn remaining_power<S: SubstrateAbsorptionSource + ?Sized>(
    source: &S,
    depth_m: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, LayeredWaveHeatingError> {
    // A passive lossy half-space ultimately absorbs its forward normal flux.
    // For a lossless exit, flux remains transmitted and is NOT heat, including
    // the exactly zero normal flux of a lossless evanescent field.
    if depth_m == 0.0 || source.original_stack().exit_permittivity_f_m().im == 0.0 {
        Ok(source.interface_power_w())
    } else {
        Ok(source.absorbed_between_w(depth_m, f64::INFINITY, cancelled)?)
    }
}

impl<'p, 'b, 'a, S: SubstrateAbsorptionSource + ?Sized>
    PreparedSubstrateWaveHeating<'p, 'b, 'a, S>
{
    /// Bind one finite depth span of the original homogeneous exit material.
    /// Every conductivity/capacity source must retain that material card and
    /// its shared fixed query context. Optical temperature remains held and
    /// explicit, as in the original layered-wave owner. The cross-sectional
    /// area must match the complete periodic cell, within roundoff only.
    pub fn try_new(
        thermal: &'b BoundLayeredThermalBody<'a>,
        source: &'p S,
        depth_start_m: f64,
        limits: SubstrateWaveHeatingLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LayeredWaveHeatingError> {
        poll(&mut cancelled)?;
        let body = thermal.prepared_body();
        let n = body.cell_count();
        if !depth_start_m.is_finite()
            || depth_start_m < 0.0
            || !body.conduction().geometry().is_planar()
            || n == 0
            || n > limits.maximum_cells
        {
            return Err(invalid(
                "substrate heat needs a finite nonnegative depth and admitted planar thermal cells",
            ));
        }
        let modes = source.mode_factor_count()?;
        let mode_cell_pairs = n
            .checked_add(2)
            .and_then(|queries| modes.checked_mul(queries))
            .ok_or(invalid("substrate cell-mode work overflows"))?;
        let per_cell = core::mem::size_of::<SubstrateCellPower>()
            + core::mem::size_of::<CellId>()
            + core::mem::size_of::<f64>()
            + core::mem::size_of::<(CellId, f64)>();
        let workspace_bytes = n
            .checked_mul(per_cell)
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<Self>()))
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<SubstrateWaveStep<'_>>()))
            .ok_or(invalid("substrate adapter storage overflows"))?;
        if mode_cell_pairs > limits.maximum_mode_cell_pairs
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(invalid("substrate heat exceeds its work or storage budget"));
        }
        let grid = source.original_grid();
        let optical_area = (Dd::from_f64(grid.width as f64)
            * Dd::from_f64(grid.spacing_m[0])
            * Dd::from_f64(grid.height as f64)
            * Dd::from_f64(grid.spacing_m[1]))
        .to_f64();
        let area = body.conduction().area_m2();
        if !optical_area.is_finite()
            || optical_area <= 0.0
            || (optical_area / area - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(invalid(
                "thermal area must equal the entire original periodic optical cell",
            ));
        }
        let original = source.original_stack().exit_material();
        for (layer, storage) in body.conduction().layers().iter().zip(body.storages()) {
            poll(&mut cancelled)?;
            if original.card_identity() != layer.material_card_identity() {
                return Err(invalid(
                    "all substrate thermal layers must retain the original exit material card",
                ));
            }
            admit_context(original, &layer.query_endpoints()[0], layer.temperature_axis())?;
            for component in storage.components() {
                poll(&mut cancelled)?;
                if original.card_identity() != component.material_card_identity() {
                    return Err(invalid(
                        "substrate capacity must retain the original exit material card",
                    ));
                }
                admit_context(
                    original,
                    &component.query_endpoints()[0],
                    component.temperature_axis(),
                )?;
            }
        }
        let ids = thermal.cell_ids();
        let mut cells = reserve(n)?;
        let mut layer_start = Dd::from_f64(depth_start_m);
        let mut left = depth_start_m;
        let mut absorbed = Dd::ZERO;
        for (layer, &count) in body.conduction().layers().iter().zip(body.layer_cell_counts()) {
            let width = Dd::from_f64(layer.thickness_m());
            for local in 0..count {
                poll(&mut cancelled)?;
                let right = (layer_start
                    + width * Dd::from_f64((local + 1) as f64) / Dd::from_f64(count as f64))
                .to_f64();
                if !right.is_finite() || right <= left {
                    return Err(invalid("substrate cell depth interval is unrepresentable"));
                }
                let watts = source.absorbed_between_w(left, right, &mut cancelled)?;
                if !watts.is_finite() || watts < 0.0 {
                    return Err(invalid("substrate cell watts are unrepresentable"));
                }
                cells.push(SubstrateCellPower {
                    cell_id: ids[cells.len()],
                    depth_m: [left, right],
                    absorbed_w: watts,
                });
                absorbed = absorbed + Dd::from_f64(watts);
                left = right;
            }
            layer_start = layer_start + width;
        }
        let entering_w = remaining_power(source, depth_start_m, &mut cancelled)?;
        let leaving_w = remaining_power(source, left, &mut cancelled)?;
        let absorbed_w = absorbed.to_f64();
        let residual_w =
            (Dd::from_f64(entering_w) - Dd::from_f64(leaving_w) - absorbed).to_f64();
        if [entering_w, leaving_w, absorbed_w].iter().any(|v| !v.is_finite() || *v < 0.0)
            || !residual_w.is_finite()
        {
            return Err(invalid("substrate power account is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            thermal,
            cells,
            power: SubstrateWavePower { entering_w, leaving_w, absorbed_w, residual_w },
            mode_cell_pairs,
            workspace_bytes,
        })
    }

    pub const fn source(&self) -> &'p S {
        self.source
    }
    pub const fn thermal_body(&self) -> &'b BoundLayeredThermalBody<'a> {
        self.thermal
    }
    pub fn cells(&self) -> &[SubstrateCellPower] {
        &self.cells
    }
    pub const fn power(&self) -> SubstrateWavePower {
        self.power
    }
    pub const fn preparation_mode_cell_pairs(&self) -> usize {
        self.mode_cell_pairs
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Advance the native thermal body under this held cycle-mean illumination.
    /// Native authored heat, boundaries and accepted histories are unchanged;
    /// exactly one additional cell-energy vector enters the original solver.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        conduction: LayerConductionConfig,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SubstrateWaveStep<'t>, LayeredWaveHeatingError> {
        advance_cells_to(self.thermal, &self.cells, self.power, initial, target_s,
            conduction, controls, workspace, cancelled)
    }
}

// Both native electromagnetic sources deposit their cached, source-bound cell
// watts through this one thermal path. It never evaluates or approximates light.
#[allow(clippy::too_many_arguments)]
fn advance_cells_to<'t>(
    body: &BoundLayeredThermalBody<'_>,
    cells: &[SubstrateCellPower],
    power: SubstrateWavePower,
    initial: &DepositionSimulation<'t>,
    target_s: f64,
    conduction: LayerConductionConfig,
    controls: LayeredThermalSolveControls,
    workspace: &mut LayeredThermalWorkspace,
    cancelled: impl FnMut() -> bool,
) -> Result<SubstrateWaveStep<'t>, LayeredWaveHeatingError> {
    let result = advance_cell_heat_to(
        body,
        cells,
        initial,
        target_s,
        conduction,
        controls,
        workspace,
        |dt, absorbed| {
            let entering_j = physical_joules(power.entering_w, dt)?;
            let leaving_j = physical_joules(power.leaving_w, dt)?;
            let absorbed_j = absorbed.to_f64();
            let residual_j =
                (Dd::from_f64(entering_j) - Dd::from_f64(leaving_j) - absorbed).to_f64();
            if !absorbed_j.is_finite() || !residual_j.is_finite() {
                return Err(invalid("substrate interval energy is unrepresentable"));
            }
            Ok(SubstrateWaveEnergy {
                entering_j,
                leaving_j,
                absorbed_j,
                residual_j,
            })
        },
        cancelled,
    )?;
    Ok(SubstrateWaveStep {
        thermal: result.thermal,
        heat_report: result.heat_report,
        optical: result.account,
        cell_energy_j: result.cell_energy_j,
    })
}

// The original thermal solve and single heat debit are independent of whether
// the electromagnetic control volume has two periodic faces or six finite
// faces. Each source owner supplies its own physical account before solving.
struct CellHeatStep<'t, A> {
    thermal: DepositionSimulation<'t>,
    heat_report: LayeredThermalSolveReport,
    cell_energy_j: Vec<(CellId, f64)>,
    account: A,
}

#[allow(clippy::too_many_arguments)]
fn advance_cell_heat_to<'t, A>(
    body: &BoundLayeredThermalBody<'_>,
    cells: &[SubstrateCellPower],
    initial: &DepositionSimulation<'t>,
    target_s: f64,
    conduction: LayerConductionConfig,
    controls: LayeredThermalSolveControls,
    workspace: &mut LayeredThermalWorkspace,
    account: impl FnOnce(f64, Dd) -> Result<A, LayeredWaveHeatingError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<CellHeatStep<'t, A>, LayeredWaveHeatingError> {
    poll(&mut cancelled)?;
    let step = LayeredThermalStep::new(body, initial, target_s, conduction, || {
        checkpoint(&mut cancelled)
    })?;
    let dt = step.duration_s();
    let mut additions = reserve(cells.len())?;
    let mut cell_energy_j = reserve(cells.len())?;
    let mut absorbed = Dd::ZERO;
    for (&id, cell) in step.cell_ids().iter().zip(cells) {
        poll(&mut cancelled)?;
        if id != cell.cell_id {
            return Err(invalid("substrate heat must retain native thermal cell order"));
        }
        let joules = physical_joules(cell.absorbed_w, dt)?;
        additions.push(joules);
        cell_energy_j.push((id, joules));
        absorbed = absorbed + Dd::from_f64(joules);
    }
    let account = account(dt, absorbed)?;
    let step = step.with_cell_energy_j(&additions, || checkpoint(&mut cancelled))?;
    let (thermal, _, heat_report) =
        step.solve(controls, workspace, || checkpoint(&mut cancelled))?;
    poll(&mut cancelled)?;
    Ok(CellHeatStep {
        thermal,
        heat_report,
        cell_energy_j,
        account,
    })
}

fn physical_joules(watts: f64, duration_s: f64) -> Result<f64, LayeredWaveHeatingError> {
    let joules = watts * duration_s;
    if !joules.is_finite() || (watts > 0.0 && joules == 0.0) {
        return Err(invalid("positive substrate energy is unrepresentable"));
    }
    Ok(joules)
}
