//! Physical, source-bound orthotropic conduction on the original thermal slabs.
//!
//! Temperature-dependent principal conductivities retain their source paths.
//! The slab's existing physical storage, SI source/force conversions, population
//! collision, wall heat accounting and acceptance transaction remain the owners.
//! Momentum diffusion and buoyancy retain the admitted background fluid model.

use super::{HeatBalance, HeatSourceLimits, MaterialSlab, ThermalLbm, ThermalLbm3};
use crate::material::electrothermal::midpoint::MAX_MIDPOINT_ITERATIONS;
use crate::thermal::sources::{Slab, ThermalSourceError, source_step_tensor_mapped};
use crate::thermal::tensor::{DiffusivityField, DiffusivityTensor};
use fs_exec::CancelGate;
use fs_material::thermal::transport::ThermalTransportError;
use fs_material::thermal::transport::orthotropic::{
    OrthotropicThermalTransportSample, PreparedOrthotropicThermalTransport,
};

/// Physical material assignments in the fixed Cartesian axes of the lattice.
/// Per-cell assignments include unused wall slots. References allow many cells
/// to reuse one original prepared owner without cloning its source receipts.
#[derive(Clone, Copy, Debug)]
pub enum ThermalMaterialField<'a, const D: usize> {
    Uniform(&'a PreparedOrthotropicThermalTransport<D>),
    PerCell(&'a [&'a PreparedOrthotropicThermalTransport<D>]),
}
impl<'a, const D: usize> ThermalMaterialField<'a, D> {
    fn at(self, cell: usize) -> &'a PreparedOrthotropicThermalTransport<D> {
        match self {
            Self::Uniform(path) => path,
            Self::PerCell(paths) => paths[cell],
        }
    }
}

/// Exact source evaluation and the discrete tensor prepared from it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialThermalTensorSample<'a, const D: usize> {
    /// Original linear fluid-cell index; wall cells are omitted.
    pub cell: usize,
    /// Original prepared source and all of its receipts, borrowed without per-cell cloning.
    pub source: &'a PreparedOrthotropicThermalTransport<D>,
    /// Source path identity, physical temperature, storage and full SI tensor.
    pub material: OrthotropicThermalTransportSample<D>,
    /// alpha_lattice = alpha_SI * dt/dx^2, using the retained physical mapping.
    pub diffusion: DiffusivityTensor<D>,
}

/// Accepted physical heat balance and both original-source evaluation epochs.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialThermalTensorStep<'a, const D: usize> {
    /// Independent physical source/storage/wall balance [J].
    pub heat: HeatBalance,
    /// Coefficients actually used for this step. The explicit entry samples
    /// its starting state; the midpoint entry samples its converged trial.
    pub applied: Vec<MaterialThermalTensorSample<'a, D>>,
    /// Endpoint support/representation checks, not an additional heat step.
    pub endpoint: Vec<MaterialThermalTensorSample<'a, D>>,
}

/// Bounded conductivity/temperature feedback over one original lattice step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorMidpointOptions {
    /// Nonlinear trials in `1..=MAX_MIDPOINT_ITERATIONS`, as in scalar feedback.
    pub max_iterations: usize,
    /// Maximum cellwise difference from the actual old/new midpoint, K.
    pub temperature_tolerance_k: f64,
    /// Maximum relative mismatch of ANY principal conductivity, in [0, 1].
    /// Each direction has its own scale, including a weak conduction axis.
    pub conductivity_relative_tolerance: f64,
    /// Under-relaxation in (0, 1]. Changes iteration, not the physical step.
    pub relaxation: f64,
}
impl TensorMidpointOptions {
    fn admit(self) -> Result<(), MaterialThermalTensorError> {
        if self.max_iterations == 0
            || self.max_iterations > MAX_MIDPOINT_ITERATIONS
            || !self.temperature_tolerance_k.is_finite()
            || self.temperature_tolerance_k <= 0.0
            || !self.conductivity_relative_tolerance.is_finite()
            || !(0.0..=1.0).contains(&self.conductivity_relative_tolerance)
            || !self.relaxation.is_finite()
            || self.relaxation <= 0.0
            || self.relaxation > 1.0
        {
            return Err(MaterialThermalTensorError::InvalidMidpointOptions);
        }
        Ok(())
    }
}

/// Original physical step and independent midpoint consistency checks.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialThermalTensorMidpointStep<'a, const D: usize> {
    pub step: MaterialThermalTensorStep<'a, D>,
    /// Attempted native steps, each starting from the same accepted state.
    pub iterations: usize,
    pub temperature_residual_k: f64,
    /// Worst directional relative mismatch at the actual old/new midpoint.
    pub conductivity_relative_residual: f64,
}

/// Every refusal leaves the complete accepted material slab unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialThermalTensorError {
    Invalid(&'static str),
    InvalidMidpointOptions,
    /// Bounded iteration failed; no partial thermal or flow step is published.
    NonConvergence {
        iterations: usize,
        temperature_residual_k: f64,
        /// None when the temperature gate prevented a constitutive recheck.
        conductivity_relative_residual: Option<f64>,
    },
    Source {
        cell: usize,
        temperature_k: f64,
        source: ThermalTransportError,
    },
    Diffusion {
        cell: usize,
        temperature_k: f64,
        source: ThermalSourceError,
    },
    Thermal(ThermalSourceError),
    Cancelled,
}
impl core::fmt::Display for MaterialThermalTensorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "material thermal tensor: {s}"),
            Self::InvalidMidpointOptions => {
                write!(f, "invalid material thermal tensor midpoint controls")
            }
            Self::NonConvergence {
                iterations,
                temperature_residual_k,
                conductivity_relative_residual,
            } => write!(
                f,
                "material thermal tensor midpoint did not converge in {iterations} trials: {temperature_residual_k} K, conductivity {conductivity_relative_residual:?}"
            ),
            Self::Source {
                cell,
                temperature_k,
                source,
            } => write!(
                f,
                "material thermal tensor at cell {cell}, {temperature_k} K: {source}"
            ),
            Self::Diffusion {
                cell,
                temperature_k,
                source,
            } => write!(
                f,
                "material thermal tensor at cell {cell}, {temperature_k} K: {source}"
            ),
            Self::Thermal(e) => write!(f, "material thermal tensor: {e}"),
            Self::Cancelled => write!(f, "material thermal tensor cancelled"),
        }
    }
}
impl std::error::Error for MaterialThermalTensorError {}
impl From<ThermalSourceError> for MaterialThermalTensorError {
    fn from(error: ThermalSourceError) -> Self {
        Self::Thermal(error)
    }
}
fn poll(gate: &CancelGate) -> Result<(), MaterialThermalTensorError> {
    if gate.is_requested() {
        Err(MaterialThermalTensorError::Cancelled)
    } else {
        Ok(())
    }
}

fn admit_storage<S, const D: usize>(
    mapping: &MaterialSlab<S>,
    field: ThermalMaterialField<'_, D>,
    cells: usize,
    plane: usize,
    gate: &CancelGate,
) -> Result<(), MaterialThermalTensorError> {
    let expected = [
        mapping.force_units.transport.density_kg_m3(),
        mapping.force_units.transport.specific_heat_j_kg_k(),
    ];
    let check = |path: &PreparedOrthotropicThermalTransport<D>| {
        if !path.constant_storage().is_some_and(|storage| {
            storage
                .iter()
                .zip(expected)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        }) {
            return Err(MaterialThermalTensorError::Invalid(
                "original rho and heat-capacity sources must equal retained storage throughout the whole temperature path",
            ));
        }
        Ok(())
    };
    match field {
        ThermalMaterialField::Uniform(path) => check(path)?,
        ThermalMaterialField::PerCell(paths) => {
            if paths.len() != cells {
                return Err(MaterialThermalTensorError::Invalid(
                    "one material assignment per cell, including unused walls, is required",
                ));
            }
            for cell in plane..cells - plane {
                if cell % 256 == 0 {
                    poll(gate)?;
                }
                check(paths[cell])?;
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn sample<'a, S: Slab<D, Q>, const D: usize, const Q: usize>(
    mapping: &MaterialSlab<()>,
    state: &S,
    field: ThermalMaterialField<'a, D>,
    cells: usize,
    plane: usize,
    limits: HeatSourceLimits,
    gate: &CancelGate,
) -> Result<Vec<MaterialThermalTensorSample<'a, D>>, MaterialThermalTensorError> {
    sample_temperatures(mapping, field, cells, plane, limits, gate, |cell| {
        let pop = &state.populations()[cell];
        if pop.iter().any(|value| !value.is_finite()) {
            return Err(MaterialThermalTensorError::Invalid(
                "nonfinite thermal populations",
            ));
        }
        mapping
            .temperature_si(pop.iter().sum())
            .map_err(|error| MaterialThermalTensorError::Thermal(error.into()))
    })
}

#[allow(clippy::too_many_arguments)]
fn sample_temperatures<'a, const D: usize>(
    mapping: &MaterialSlab<()>,
    field: ThermalMaterialField<'a, D>,
    cells: usize,
    plane: usize,
    limits: HeatSourceLimits,
    gate: &CancelGate,
    mut temperature: impl FnMut(usize) -> Result<f64, MaterialThermalTensorError>,
) -> Result<Vec<MaterialThermalTensorSample<'a, D>>, MaterialThermalTensorError> {
    let units = mapping.force_units;
    let scale = (units.timestep_s / units.cell_size_m) / units.cell_size_m;
    if !scale.is_finite() || scale <= 0.0 {
        return Err(MaterialThermalTensorError::Invalid(
            "diffusivity scale is unrepresentable",
        ));
    }
    let mut samples = Vec::with_capacity(cells - 2 * plane);
    // Reuse a factorization when adjacent cells have identical evaluated
    // tensors. In particular uniform constant laws require one preparation,
    // even when the temperature field varies. Source evaluations stay local.
    let mut previous: Option<([[f64; D]; D], DiffusivityTensor<D>)> = None;
    for cell in plane..cells - plane {
        if cell % 256 == 0 {
            poll(gate)?;
        }
        let temperature_k = temperature(cell)?;
        if !temperature_k.is_finite() || temperature_k < limits.temperature_range_k[0]
            || temperature_k > limits.temperature_range_k[1]
        {
            return Err(MaterialThermalTensorError::Invalid(
                "physical temperature is outside the admitted interval",
            ));
        }
        let material = field
            .at(cell)
            .at_temperature(temperature_k)
            .map_err(|source| MaterialThermalTensorError::Source {
                cell,
                temperature_k,
                source,
            })?;
        let physical = material.diffusivity_m2_s();
        let mut lattice = [[0.0; D]; D];
        for row in 0..D {
            for column in 0..D {
                lattice[row][column] = physical[row][column] * scale;
                if !lattice[row][column].is_finite()
                    || (physical[row][column] != 0.0 && lattice[row][column] == 0.0)
                {
                    return Err(MaterialThermalTensorError::Invalid(
                        "physical tensor is unrepresentable in lattice units",
                    ));
                }
            }
        }
        let diffusion = match previous {
            Some((tensor, prepared)) if tensor == lattice => prepared,
            _ => DiffusivityTensor::from_lattice(lattice).map_err(|source| {
                MaterialThermalTensorError::Diffusion {
                    cell,
                    temperature_k,
                    source,
                }
            })?,
        };
        previous = Some((lattice, diffusion));
        samples.push(MaterialThermalTensorSample {
            cell,
            source: field.at(cell),
            material,
            diffusion,
        });
    }
    poll(gate)?;
    Ok(samples)
}

#[allow(clippy::too_many_arguments)]
fn advance<'a, S: Slab<D, Q>, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>,
    power: &[f64],
    force: &[[f64; D]],
    depth_m: f64,
    field: ThermalMaterialField<'a, D>,
    limits: HeatSourceLimits,
    midpoint: Option<TensorMidpointOptions>,
    gate: &CancelGate,
) -> Result<MaterialThermalTensorMidpointStep<'a, D>, MaterialThermalTensorError> {
    poll(gate)?;
    if let Some(options) = midpoint {
        options.admit()?;
    }
    let (cells, plane) = slab.lattice.layout()?;
    if cells > limits.max_cells
        || slab.lattice.walls() != slab.wall_coordinates
        || !slab.lattice.flow_valid()
    {
        return Err(MaterialThermalTensorError::Invalid(
            "invalid physical slab, wall mapping or cell budget",
        ));
    }
    // Copy only the retained physical mapping. A unit placeholder stores no
    // populations and lets temperature_si/heat_inputs run unchanged, while
    // the native transaction exclusively borrows the actual solver.
    let mapping = MaterialSlab {
        lattice: (),
        top_temperature_k: slab.top_temperature_k,
        temperature_difference_k: slab.temperature_difference_k,
        wall_coordinates: slab.wall_coordinates,
        force_units: slab.force_units,
    };
    admit_storage(&mapping, field, cells, plane, gate)?;
    let input = mapping.heat_inputs(power, force, cells, depth_m, limits)?;
    if power[..plane]
        .iter()
        .chain(&power[cells - plane..])
        .any(|x| *x != 0.0)
    {
        return Err(MaterialThermalTensorError::Invalid(
            "wall power must be zero",
        ));
    }
    let mut applied =
        sample::<S, D, Q>(&mapping, &slab.lattice, field, cells, plane, limits, gate)?;
    let old: Vec<_> = if midpoint.is_some() {
        applied
            .iter()
            .map(|sample| sample.material.temperature_k())
            .collect()
    } else {
        Vec::new()
    };
    let mut evaluation = old.clone();
    // Both native slab layouts have at least two fluid planes. Wall tensor
    // slots reuse an admitted fluid value; wall material laws are never queried.
    let mut tensors = vec![applied[0].diffusion; cells];
    let mut iterations = 0;
    loop {
        iterations += 1;
        for value in &applied {
            tensors[value.cell] = value.diffusion;
        }
        poll(gate)?;
        // The native source owner snapshots once per trial and adopts only
        // after this callback succeeds. A midpoint retry uses that SAME
        // refusal boundary, so no extra grid snapshot or accepted prefix is
        // needed. Histories, force frames and physical heat all restart together.
        let mut endpoint_refusal = None;
        let result = source_step_tensor_mapped::<S, D, Q, _>(
            &mut slab.lattice,
            &input.heat,
            &input.force,
            input.limits,
            DiffusivityField::PerCell(&tensors),
            |candidate, balance| {
                let heat = input.balance(balance)?;
                let checked = (|| {
                    let endpoint =
                        sample::<S, D, Q>(&mapping, candidate, field, cells, plane, limits, gate)?;
                    let mut temperature_residual_k = 0.0_f64;
                    let mut conductivity_relative_residual = 0.0_f64;
                    if let Some(options) = midpoint {
                        for (index, sample) in endpoint.iter().enumerate() {
                            if index % 256 == 0 {
                                poll(gate)?;
                            }
                            let actual = 0.5 * old[index] + 0.5 * sample.material.temperature_k();
                            temperature_residual_k =
                                temperature_residual_k.max((actual - evaluation[index]).abs());
                        }
                        let coordinate_ok =
                            temperature_residual_k <= options.temperature_tolerance_k;
                        if coordinate_ok {
                            conductivity_relative_residual =
                                conductivity_residual(&applied, &endpoint, &old, gate)?;
                        }
                        if !coordinate_ok
                            || conductivity_relative_residual
                                > options.conductivity_relative_tolerance
                        {
                            for (index, sample) in endpoint.iter().enumerate() {
                                if index % 256 == 0 {
                                    poll(gate)?;
                                }
                                let actual =
                                    0.5 * old[index] + 0.5 * sample.material.temperature_k();
                                evaluation[index] = options
                                    .relaxation
                                    .mul_add(actual - evaluation[index], evaluation[index]);
                            }
                            return Err(MaterialThermalTensorError::NonConvergence {
                                iterations,
                                temperature_residual_k,
                                conductivity_relative_residual: coordinate_ok
                                    .then_some(conductivity_relative_residual),
                            });
                        }
                    }
                    poll(gate)?;
                    Ok((
                        endpoint,
                        temperature_residual_k,
                        conductivity_relative_residual,
                    ))
                })()
                .map_err(|error| {
                    endpoint_refusal = Some(error);
                    ThermalSourceError::Invalid("material tensor endpoint refused")
                })?;
                Ok(MaterialThermalTensorMidpointStep {
                    step: MaterialThermalTensorStep {
                        heat,
                        applied,
                        endpoint: checked.0,
                    },
                    iterations,
                    temperature_residual_k: checked.1,
                    conductivity_relative_residual: checked.2,
                })
            },
        );
        match result {
            Ok(step) => return Ok(step),
            Err(error) => {
                let error = endpoint_refusal.unwrap_or_else(|| error.into());
                if !matches!(error, MaterialThermalTensorError::NonConvergence { .. })
                    || iterations >= midpoint.map_or(1, |options| options.max_iterations)
                {
                    return Err(error);
                }
                applied =
                    sample_temperatures(&mapping, field, cells, plane, limits, gate, |cell| {
                        Ok(evaluation[cell - plane])
                    })?;
            }
        }
    }
}

// Compare the original positive principal laws, not a matrix norm dominated
// by the strongest direction or a cancellation-prone rotation back from K.
// Fixed rho*Cp makes this the same relative check for physical diffusivity.
fn conductivity_residual<const D: usize>(
    applied: &[MaterialThermalTensorSample<'_, D>],
    endpoint: &[MaterialThermalTensorSample<'_, D>],
    old: &[f64],
    gate: &CancelGate,
) -> Result<f64, MaterialThermalTensorError> {
    let mut residual = 0.0_f64;
    for (index, (used, end)) in applied.iter().zip(endpoint).enumerate() {
        if index % 256 == 0 {
            poll(gate)?;
        }
        let actual = 0.5 * old[index] + 0.5 * end.material.temperature_k();
        for path in used.source.principal_paths() {
            let conductivity = |temperature_k| {
                path.at_temperature(temperature_k)
                    .map(|value| value.conductivity_w_m_k())
                    .map_err(|source| MaterialThermalTensorError::Source {
                        cell: used.cell,
                        temperature_k,
                        source,
                    })
            };
            let a = conductivity(used.material.temperature_k())?;
            let b = conductivity(actual)?;
            residual = residual.max((a - b).abs() / a.max(b));
        }
    }
    Ok(residual)
}

impl MaterialSlab<ThermalLbm> {
    /// Advance physical heat and flow using source-bound, oriented K(T).
    ///
    /// Power is W/m3, external force is N/m3, and explicit depth is metres.
    /// Every principal source must share its original density/Cp sources and
    /// context. Both storage laws must be constant over their WHOLE admitted
    /// paths and match the slab's retained storage exactly. Variable storage,
    /// enthalpy/latent heat and material contact resistance require other owners.
    ///
    /// Principal axes are expressed in the square lattice frame. Coefficients
    /// are held at the initial cell temperatures for one interval; call this
    /// method each step for feedback. Endpoint source and tensor checks precede
    /// adoption. This is explicit feedback, not a nonlinear implicit solution.
    /// Viscosity, even thermal moments and buoyancy keep their existing setup;
    /// scalar background diffusivity/Nusselt diagnostics do not describe K.
    ///
    /// Cancellation is checked around bounded source loops and whole native
    /// kernel steps. Source, SI balance, endpoint and cancellation failures
    /// preserve the original populations, force state and advection history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_tensor_diffusion<'a>(
        &mut self,
        power_w_m3: &[f64],
        force_n_m3: &[[f64; 2]],
        depth_m: f64,
        materials: ThermalMaterialField<'a, 2>,
        limits: HeatSourceLimits,
        gate: &CancelGate,
    ) -> Result<MaterialThermalTensorStep<'a, 2>, MaterialThermalTensorError> {
        advance::<ThermalLbm, 2, 5>(
            self, power_w_m3, force_n_m3, depth_m, materials, limits, None, gate,
        )
        .map(|result| result.step)
    }

    /// Couple full source-bound K(T) to the actual old/new temperature midpoint.
    ///
    /// Solves `state_next = native_step(state_old, K((T_old + T_next)/2))`.
    /// Every trial uses the original source kicks, tensor collision, streaming,
    /// wall exchange and flow step. Both temperature consistency and freshly
    /// checked principal conductivities must pass before the only publication.
    /// Source/endpoint/energy refusal, cancellation or exhausted work preserves
    /// all original populations, physical forces and advection history.
    ///
    /// Power and external force are prescribed over this same interval. The
    /// original frozen-storage and material-domain restrictions still apply.
    /// This removes conductivity lag from the declared fixed-point equation;
    /// it does not establish second-order accuracy or unconditional stability
    /// of the complete split flow/source solver. Refine time and space separately.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_tensor_diffusion_midpoint<'a>(
        &mut self,
        power_w_m3: &[f64],
        force_n_m3: &[[f64; 2]],
        depth_m: f64,
        materials: ThermalMaterialField<'a, 2>,
        limits: HeatSourceLimits,
        options: TensorMidpointOptions,
        gate: &CancelGate,
    ) -> Result<MaterialThermalTensorMidpointStep<'a, 2>, MaterialThermalTensorError> {
        advance::<ThermalLbm, 2, 5>(
            self,
            power_w_m3,
            force_n_m3,
            depth_m,
            materials,
            limits,
            Some(options),
            gate,
        )
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// D3Q19/D3Q7 counterpart with physical cell volume dx^3 and full 3D
    /// principal-to-lattice rotation. All source and acceptance rules are shared.
    pub fn step_with_material_tensor_diffusion<'a>(
        &mut self,
        power_w_m3: &[f64],
        force_n_m3: &[[f64; 3]],
        materials: ThermalMaterialField<'a, 3>,
        limits: HeatSourceLimits,
        gate: &CancelGate,
    ) -> Result<MaterialThermalTensorStep<'a, 3>, MaterialThermalTensorError> {
        let depth_m = self.force_units.cell_size_m;
        advance::<ThermalLbm3, 3, 7>(
            self, power_w_m3, force_n_m3, depth_m, materials, limits, None, gate,
        )
        .map(|result| result.step)
    }

    /// D3Q19/D3Q7 midpoint counterpart, using the original physical dx^3 cells.
    /// The same bounded nonlinear iteration, directional checks and transaction
    /// serve both dimensions; no directional average substitutes for the tensor.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_tensor_diffusion_midpoint<'a>(
        &mut self,
        power_w_m3: &[f64],
        force_n_m3: &[[f64; 3]],
        materials: ThermalMaterialField<'a, 3>,
        limits: HeatSourceLimits,
        options: TensorMidpointOptions,
        gate: &CancelGate,
    ) -> Result<MaterialThermalTensorMidpointStep<'a, 3>, MaterialThermalTensorError> {
        let depth_m = self.force_units.cell_size_m;
        advance::<ThermalLbm3, 3, 7>(
            self,
            power_w_m3,
            force_n_m3,
            depth_m,
            materials,
            limits,
            Some(options),
            gate,
        )
    }
}
