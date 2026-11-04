//! Physical, source-bound orthotropic conduction on the original thermal slabs.
//!
//! Temperature-dependent principal conductivities retain their source paths.
//! The slab's existing physical storage, SI source/force conversions, population
//! collision, wall heat accounting and acceptance transaction remain the owners.
//! Momentum diffusion and buoyancy retain the admitted background fluid model.

use super::{HeatBalance, HeatSourceLimits, MaterialSlab, ThermalLbm, ThermalLbm3};
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
    /// Coefficients used for this step, held at its accepted starting state.
    pub applied: Vec<MaterialThermalTensorSample<'a, D>>,
    /// Endpoint support/representation checks, not an additional heat step.
    pub endpoint: Vec<MaterialThermalTensorSample<'a, D>>,
}

/// Every refusal leaves the complete accepted material slab unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialThermalTensorError {
    Invalid(&'static str),
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
        let pop = &state.populations()[cell];
        if pop.iter().any(|value| !value.is_finite()) {
            return Err(MaterialThermalTensorError::Invalid(
                "nonfinite thermal populations",
            ));
        }
        let temperature_k = mapping
            .temperature_si(pop.iter().sum())
            .map_err(ThermalSourceError::from)?;
        if temperature_k < limits.temperature_range_k[0]
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
    gate: &CancelGate,
) -> Result<MaterialThermalTensorStep<'a, D>, MaterialThermalTensorError> {
    poll(gate)?;
    let (cells, plane) = slab.lattice.layout()?;
    if cells > limits.max_cells || slab.lattice.walls() != (1.0, 0.0) || !slab.lattice.flow_valid()
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
    let applied = sample::<S, D, Q>(&mapping, &slab.lattice, field, cells, plane, limits, gate)?;
    // Both native slab layouts have at least two fluid planes. Wall tensor
    // slots reuse an admitted fluid value; wall material laws are never queried.
    let mut tensors = vec![applied[0].diffusion; cells];
    for value in &applied {
        tensors[value.cell] = value.diffusion;
    }
    poll(gate)?;
    // The native source owner snapshots ONCE and adopts only after this
    // callback succeeds. Keep the richer material/cancellation refusal intact
    // across its thermal-only callback error type; no outer grid snapshot is needed.
    let mut endpoint_refusal = None;
    let result = source_step_tensor_mapped::<S, D, Q, _>(
        &mut slab.lattice,
        &input.heat,
        &input.force,
        input.limits,
        DiffusivityField::PerCell(&tensors),
        |candidate, balance| {
            let heat = input.balance(balance)?;
            let endpoint =
                sample::<S, D, Q>(&mapping, candidate, field, cells, plane, limits, gate)
                    .and_then(|endpoint| {
                        poll(gate)?;
                        Ok(endpoint)
                    })
                    .map_err(|error| {
                        endpoint_refusal = Some(error);
                        ThermalSourceError::Invalid("material tensor endpoint refused")
                    })?;
            Ok((heat, endpoint))
        },
    );
    let (heat, endpoint) =
        result.map_err(|error| endpoint_refusal.unwrap_or_else(|| error.into()))?;
    Ok(MaterialThermalTensorStep {
        heat,
        applied,
        endpoint,
    })
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
            self, power_w_m3, force_n_m3, depth_m, materials, limits, gate,
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
            self, power_w_m3, force_n_m3, depth_m, materials, limits, gate,
        )
    }
}
