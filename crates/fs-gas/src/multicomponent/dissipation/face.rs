//! Shared original-source molecular faces for axial, planar and spatial fields.
use super::evolve::{close_simplex, fractions, heat_capacity_density};
use super::law::validate_binary;
use super::*;
use crate::multicomponent::{MAX_NORMAL_COMPONENTS, MAX_SPECIES, NormalCell, NormalPrimitive, sum};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

mod spatial;
pub use spatial::SpatialTransport;

/// Original constitutive residuals and a conservative engineering diffusion scale.
/// The scale includes composition/thermal coupling and the complete Newtonian
/// viscous tensor, but is not a spectral or invariant-domain certificate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransportFaceReport {
    /// Maximum face diffusion scale before spatial inverse lengths [m²/s].
    pub diffusivity_m2_s: f64,
    /// Largest checked original Maxwell--Stefan row residual [1/m].
    pub maximum_gradient_residual_per_m: f64,
    /// Signed sum of all diffusive species mass fluxes [kg/(m² s)].
    pub mass_frame_residual_kg_m2_s: f64,
}

/// Prepared molecular faces on the exact original evolving species bank.
///
/// The same private kernel serves the original axial and network owners.
/// Planar stress uses the three-dimensional Newtonian trace convention with
/// zero out-of-plane velocity, including tangential velocity derivatives.
/// Species enthalpy, Fourier heat and stress work share one conservative energy
/// flux. No separate reaction or viscous heat source is added.
///
/// Face scratch is reserved once; successful kernels allocate no new scratch.
/// The supplied law keeps its own allocation/validity contract: in particular,
/// a material-card provider may perform allocating source queries. Soret/Dufour
/// effects, no-slip walls and nonideal gas transport require other closures.
pub struct PlanarTransport<L: TransportLaw> {
    gas: MixtureGas,
    law: L,
    diffusion: MaxwellStefanV1,
    mass_ratio: f64,
    binary: Vec<f64>,
    workspace: MaxwellStefanWorkspaceV1,
}

impl<L: TransportLaw> PlanarTransport<L> {
    /// Bind identical original caloric species and one molecular transport law.
    pub fn new(
        gas: MixtureGas,
        law: L,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if law.source() != gas.source() {
            return Err(DissipationError::Invalid(
                "planar molecular law and caloric species differ",
            ));
        }
        let (diffusion, mass_ratio) = prepare_diffusion(&gas)?;
        let binary = allocated(gas.species_count() * gas.species_count(), 0.0)?;
        let workspace = diffusion.workspace()?;
        poll(&mut cancelled)?;
        Ok(Self {
            gas,
            law,
            diffusion,
            mass_ratio,
            binary,
            workspace,
        })
    }

    /// Exact original caloric source and admitted temperature interval.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Original law, including its source identities and physical domain.
    pub const fn law(&self) -> &L {
        &self.law
    }
    /// Prepared face/Maxwell--Stefan numerical heap capacity. Immutable gas and
    /// species storage, the caller-supplied law's storage/evaluation allocations,
    /// fixed-size owners and bounded stack scratch are outside this budget.
    pub fn workspace_bytes(&self) -> usize {
        self.binary.capacity() * core::mem::size_of::<f64>() + self.workspace.workspace_bytes()
    }

    /// Conservative normal flux `[j_i..., -tau_nn, -tau_tn, E_flux]`.
    /// The center-to-center distance is positive. Tangential derivatives are
    /// `[d_t u_n, d_t u_t]` in the same orthonormal face frame. Output is scratch
    /// until success; neither input trace is modified.
    #[allow(clippy::too_many_arguments)]
    pub fn face(
        &mut self,
        left: NormalCell<'_>,
        right: NormalCell<'_>,
        distance: f64,
        tangential_velocity_gradient: [f64; 2],
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let left = conservative_trace(&self.gas, left, &mut cancelled)?;
        let right = conservative_trace(&self.gas, right, &mut cancelled)?;
        interior::<true, L>(
            &self.gas,
            &mut self.law,
            &self.diffusion,
            &mut self.workspace,
            left,
            right,
            distance,
            distance,
            tangential_velocity_gradient,
            self.mass_ratio,
            &mut self.binary,
            out,
            controls,
            &mut cancelled,
        )
    }

    /// Re-admit actual reconstructed primitives through the original caloric
    /// owner without repeating its bounded energy inversion. Spatial stages
    /// that already computed velocity/T/p can use this path on every face.
    #[allow(clippy::too_many_arguments)]
    pub fn face_from_primitives(
        &mut self,
        left: NormalPrimitive<'_>,
        right: NormalPrimitive<'_>,
        distance: f64,
        tangential_velocity_gradient: [f64; 2],
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let left = primitive_trace(&self.gas, left, &mut cancelled)?;
        let right = primitive_trace(&self.gas, right, &mut cancelled)?;
        interior::<true, L>(
            &self.gas,
            &mut self.law,
            &self.diffusion,
            &mut self.workspace,
            left,
            right,
            distance,
            distance,
            tangential_velocity_gradient,
            self.mass_ratio,
            &mut self.binary,
            out,
            controls,
            &mut cancelled,
        )
    }

    /// Held exterior state at a physical reservoir face. `boundary_on_left`
    /// selects which of the two admitted traces is that prescribed state.
    /// Species diffusion, traction and Fourier coefficients retain the shared
    /// two-trace midpoint approximation. Species energy uses the selected
    /// bath's original referenced enthalpies and stress work uses its velocity:
    /// a stationary bath never receives fictitious mechanical work from a
    /// cell-center midpoint velocity. Distance is cell-center to boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn reservoir_face_from_primitives(
        &mut self,
        left: NormalPrimitive<'_>,
        right: NormalPrimitive<'_>,
        distance: f64,
        tangential_velocity_gradient: [f64; 2],
        boundary_on_left: bool,
        controls: MaxwellStefanSolveV1,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let left = primitive_trace(&self.gas, left, &mut cancelled)?;
        let right = primitive_trace(&self.gas, right, &mut cancelled)?;
        let left_velocity = [left.transverse_velocity];
        let right_velocity = [right.transverse_velocity];
        let left = left.with_transverse(left_velocity);
        let right = right.with_transverse(right_velocity);
        let effort = if boundary_on_left {
            left.effort()
        } else {
            right.effort()
        };
        interior_components(
            &self.gas,
            &mut self.law,
            &self.diffusion,
            &mut self.workspace,
            left,
            right,
            distance,
            distance,
            [[
                tangential_velocity_gradient[0],
                tangential_velocity_gradient[1],
                0.0,
            ]],
            Some(effort),
            self.mass_ratio,
            &mut self.binary,
            out,
            controls,
            &mut cancelled,
        )
    }

    /// Stationary adiabatic impermeable slip wall. Its normal traction includes
    /// `d_n u_n` and `d_t u_t`; species, shear traction and total-energy flux are
    /// exactly zero. Derivatives use the globally oriented face axes, including
    /// the appropriate lower/upper sign of the normal wall velocity difference.
    pub fn wall(
        &mut self,
        cell: NormalCell<'_>,
        normal_velocity_gradient: f64,
        transverse_tangential_gradient: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let trace = conservative_trace(&self.gas, cell, &mut cancelled)?;
        wall::<true, L>(
            &self.gas,
            &mut self.law,
            trace,
            normal_velocity_gradient,
            transverse_tangential_gradient,
            &mut self.binary,
            out,
            &mut cancelled,
        )
    }

    /// Same slip-wall law from source-admitted reconstructed primitives.
    pub fn wall_from_primitive(
        &mut self,
        cell: NormalPrimitive<'_>,
        normal_velocity_gradient: f64,
        transverse_tangential_gradient: f64,
        out: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let trace = primitive_trace(&self.gas, cell, &mut cancelled)?;
        wall::<true, L>(
            &self.gas,
            &mut self.law,
            trace,
            normal_velocity_gradient,
            transverse_tangential_gradient,
            &mut self.binary,
            out,
            &mut cancelled,
        )
    }
}

pub(super) fn prepare_diffusion(
    gas: &MixtureGas,
) -> Result<(MaxwellStefanV1, f64), DissipationError> {
    let species = gas
        .source()
        .models()
        .iter()
        .map(|model| (model.species().clone(), model.molar_mass()))
        .collect();
    let diffusion = MaxwellStefanV1::new(species)?;
    let minimum = gas
        .source()
        .models()
        .iter()
        .map(|m| m.molar_mass().value())
        .fold(f64::INFINITY, f64::min);
    let maximum = gas
        .source()
        .models()
        .iter()
        .map(|m| m.molar_mass().value())
        .fold(0.0, f64::max);
    let mass_ratio = maximum / minimum;
    if !positive(mass_ratio) {
        return Err(DissipationError::Invalid(
            "molecular mass ratio is unrepresentable",
        ));
    }
    Ok((diffusion, mass_ratio))
}

pub(super) struct Trace<'a> {
    pub densities: &'a [f64],
    pub primitive: Primitive,
    pub transverse_velocity: f64,
    pub cv_density: f64,
}

struct VectorTrace<'a, const TRANSVERSE: usize> {
    densities: &'a [f64],
    primitive: Primitive,
    transverse_velocity: [f64; TRANSVERSE],
    cv_density: f64,
}

/// Original source-admitted energy effort at a physical boundary.
#[derive(Clone, Copy)]
struct EnergyEffort<const TRANSVERSE: usize> {
    temperature: f64,
    normal_velocity: f64,
    transverse_velocity: [f64; TRANSVERSE],
}
impl<const TRANSVERSE: usize> VectorTrace<'_, TRANSVERSE> {
    fn effort(&self) -> EnergyEffort<TRANSVERSE> {
        EnergyEffort {
            temperature: self.primitive.temperature,
            normal_velocity: self.primitive.velocity,
            transverse_velocity: self.transverse_velocity,
        }
    }
}

impl<'a> Trace<'a> {
    fn with_transverse<const TRANSVERSE: usize>(
        self,
        transverse_velocity: [f64; TRANSVERSE],
    ) -> VectorTrace<'a, TRANSVERSE> {
        VectorTrace {
            densities: self.densities,
            primitive: self.primitive,
            transverse_velocity,
            cv_density: self.cv_density,
        }
    }
}

fn conservative_trace<'a>(
    gas: &MixtureGas,
    cell: NormalCell<'a>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trace<'a>, DissipationError> {
    let primitive = gas.normal_primitive(cell, cancelled)?;
    let cv_density = heat_capacity_density(gas, cell.partial_densities, primitive)?;
    Ok(Trace {
        densities: cell.partial_densities,
        primitive,
        transverse_velocity: cell.transverse_momentum / primitive.density,
        cv_density,
    })
}

fn primitive_trace<'a>(
    gas: &MixtureGas,
    cell: NormalPrimitive<'a>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Trace<'a>, DissipationError> {
    poll(cancelled)?;
    let mut q = [0.0; MAX_NORMAL_COMPONENTS];
    let primitive = gas.normal_conserved_into(cell, &mut q[..gas.species_count() + 3])?;
    let cv_density = heat_capacity_density(gas, cell.partial_densities, primitive)?;
    Ok(Trace {
        densities: cell.partial_densities,
        primitive,
        transverse_velocity: cell.transverse_velocity,
        cv_density,
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn interior<const PLANAR: bool, L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    diffusion: &MaxwellStefanV1,
    workspace: &mut MaxwellStefanWorkspaceV1,
    left: Trace<'_>,
    right: Trace<'_>,
    wl: f64,
    wr: f64,
    tangential_gradient: [f64; 2],
    mass_ratio: f64,
    binary: &mut [f64],
    out: &mut [f64],
    controls: MaxwellStefanSolveV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TransportFaceReport, DissipationError> {
    if PLANAR {
        let transverse_left = [left.transverse_velocity];
        let transverse_right = [right.transverse_velocity];
        interior_components(
            gas,
            law,
            diffusion,
            workspace,
            left.with_transverse(transverse_left),
            right.with_transverse(transverse_right),
            wl,
            wr,
            [[tangential_gradient[0], tangential_gradient[1], 0.0]],
            None,
            mass_ratio,
            binary,
            out,
            controls,
            cancelled,
        )
    } else {
        interior_components(
            gas,
            law,
            diffusion,
            workspace,
            left.with_transverse([]),
            right.with_transverse([]),
            wl,
            wr,
            [],
            None,
            mass_ratio,
            binary,
            out,
            controls,
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn interior_components<const TRANSVERSE: usize, L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    diffusion: &MaxwellStefanV1,
    workspace: &mut MaxwellStefanWorkspaceV1,
    left: VectorTrace<'_, TRANSVERSE>,
    right: VectorTrace<'_, TRANSVERSE>,
    wl: f64,
    wr: f64,
    tangential_gradient: [[f64; 3]; TRANSVERSE],
    boundary_effort: Option<EnergyEffort<TRANSVERSE>>,
    mass_ratio: f64,
    binary: &mut [f64],
    out: &mut [f64],
    controls: MaxwellStefanSolveV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TransportFaceReport, DissipationError> {
    poll(cancelled)?;
    let n = gas.species_count();
    if out.len() != n + 2 + TRANSVERSE
        || !positive(wl)
        || !positive(wr)
        || !tangential_gradient.iter().flatten().all(|v| v.is_finite())
    {
        return Err(DissipationError::Invalid(
            "molecular face dimensions or gradients",
        ));
    }
    let pl = left.primitive;
    let pr = right.primitive;
    let distance = 0.5 * wl + 0.5 * wr;
    let right_weight = (0.5 * wl) / distance;
    let left_weight = 1.0 - right_weight;
    let t = left_weight * pl.temperature + right_weight * pr.temperature;
    let p = left_weight * pl.pressure + right_weight * pr.pressure;
    let velocity = boundary_effort.map_or_else(
        || left_weight * pl.velocity + right_weight * pr.velocity,
        |effort| effort.normal_velocity,
    );
    let enthalpy_temperature = boundary_effort.map_or(t, |effort| effort.temperature);
    let grad_t = (pr.temperature - pl.temperature) / distance;
    let grad_u = (pr.velocity - pl.velocity) / distance;
    let mut xl = [0.0; MAX_SPECIES];
    let mut xr = [0.0; MAX_SPECIES];
    let mut x = [0.0; MAX_SPECIES];
    let mut gradient = [0.0; MAX_SPECIES];
    fractions(gas, left.densities, &mut xl[..n])?;
    fractions(gas, right.densities, &mut xr[..n])?;
    for k in 0..n {
        x[k] = left_weight * xl[k] + right_weight * xr[k];
        gradient[k] = (xr[k] - xl[k]) / distance;
    }
    let pivot = (0..n)
        .max_by(|&a, &b| x[a].total_cmp(&x[b]))
        .expect("nonempty species bank");
    close_simplex(&mut x[..n], pivot)?;
    gradient[pivot] = -sum((0..n).filter(|&k| k != pivot).map(|k| gradient[k]));
    let ratio = pr.pressure / pl.pressure;
    let log_pressure = if positive(ratio) {
        fs_math::det::ln(ratio)
    } else {
        fs_math::det::ln(pr.pressure) - fs_math::det::ln(pl.pressure)
    };
    let coefficients = law.evaluate(
        Temperature::new(t),
        Pressure::new(p),
        &x[..n],
        binary,
        &mut || cancelled(),
    )?;
    coefficients.validate()?;
    let maximum_diffusion = validate_binary(binary, n)?;
    let evaluated = diffusion.evaluate(
        &x[..n],
        &gradient[..n],
        Temperature::new(t),
        Pressure::new(p),
        log_pressure / distance,
        binary,
        controls,
        workspace,
        &mut *cancelled,
    )?;
    let mut report = TransportFaceReport {
        diffusivity_m2_s: 0.0,
        maximum_gradient_residual_per_m: evaluated.maximum_gradient_residual_per_m(),
        mass_frame_residual_kg_m2_s: evaluated.mass_frame_residual_kg_m2_s(),
    };
    out[..n].copy_from_slice(evaluated.mass_fluxes());
    let longitudinal = coefficients.longitudinal_viscosity()?;
    let mu = coefficients.dynamic_viscosity.value();
    let cross = coefficients.bulk_viscosity.value() - (2.0 / 3.0) * mu;
    let mut stress = longitudinal * grad_u;
    let mut shear_work = 0.0;
    for (axis, gradient) in tangential_gradient.iter().enumerate() {
        stress += cross * gradient[axis + 1];
        let transverse_gradient =
            (right.transverse_velocity[axis] - left.transverse_velocity[axis]) / distance;
        let shear = mu * (transverse_gradient + gradient[0]);
        out[n + 1 + axis] = -shear;
        let work_velocity = boundary_effort.map_or_else(
            || {
                left_weight * left.transverse_velocity[axis]
                    + right_weight * right.transverse_velocity[axis]
            },
            |effort| effort.transverse_velocity[axis],
        );
        shear_work += work_velocity * shear;
    }
    out[n] = -stress;
    let mut enthalpy_flux = [0.0; MAX_SPECIES];
    for k in 0..n {
        if out[k] == 0.0 {
            continue;
        }
        let model = &gas.source().models()[k];
        let h_referenced = (model
            .molar_internal_energy_change(
                Temperature::new(gas.domain().minimum_temperature),
                Temperature::new(enthalpy_temperature),
            )
            .map_err(FlowError::from)?
            .value()
            + R * enthalpy_temperature)
            / model.molar_mass().value();
        enthalpy_flux[k] = h_referenced * out[k];
    }
    let energy = sum(enthalpy_flux[..n].iter().copied())
        - coefficients.thermal_conductivity.value() * grad_t
        - velocity * stress;
    out[n + 1 + TRANSVERSE] = if TRANSVERSE > 0 {
        energy - shear_work
    } else {
        energy
    };
    if !out.iter().all(|value| value.is_finite()) {
        return Err(DissipationError::Invalid("non-finite molecular face flux"));
    }
    let gamma = 1.0
        + (pl.pressure / (pl.temperature * left.cv_density))
            .max(pr.pressure / (pr.temperature * right.cv_density));
    let viscosity_scale = viscosity_scale::<TRANSVERSE>(longitudinal, cross, mu);
    report.diffusivity_m2_s = (maximum_diffusion * mass_ratio * gamma)
        .max(viscosity_scale / pl.density.min(pr.density))
        .max(coefficients.thermal_conductivity.value() / left.cv_density.min(right.cv_density));
    if !report.diffusivity_m2_s.is_finite() || report.diffusivity_m2_s < 0.0 {
        return Err(DissipationError::Invalid("molecular face scale overflow"));
    }
    poll(cancelled)?;
    Ok(report)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn wall<const PLANAR: bool, L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    trace: Trace<'_>,
    normal_gradient: f64,
    transverse_gradient: f64,
    binary: &mut [f64],
    out: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TransportFaceReport, DissipationError> {
    if PLANAR {
        let transverse_velocity = [trace.transverse_velocity];
        wall_components(
            gas,
            law,
            trace.with_transverse(transverse_velocity),
            normal_gradient,
            [transverse_gradient],
            binary,
            out,
            cancelled,
        )
    } else {
        wall_components(
            gas,
            law,
            trace.with_transverse([]),
            normal_gradient,
            [],
            binary,
            out,
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn wall_components<const TRANSVERSE: usize, L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    trace: VectorTrace<'_, TRANSVERSE>,
    normal_gradient: f64,
    transverse_gradients: [f64; TRANSVERSE],
    binary: &mut [f64],
    out: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TransportFaceReport, DissipationError> {
    poll(cancelled)?;
    let n = gas.species_count();
    if out.len() != n + 2 + TRANSVERSE
        || !normal_gradient.is_finite()
        || !transverse_gradients.iter().all(|v| v.is_finite())
    {
        return Err(DissipationError::Invalid(
            "molecular wall dimensions or gradients",
        ));
    }
    let p = trace.primitive;
    let mut x = [0.0; MAX_SPECIES];
    fractions(gas, trace.densities, &mut x[..n])?;
    let c = law.evaluate(
        Temperature::new(p.temperature),
        Pressure::new(p.pressure),
        &x[..n],
        binary,
        &mut || cancelled(),
    )?;
    c.validate()?;
    validate_binary(binary, n)?;
    let longitudinal = c.longitudinal_viscosity()?;
    let mu = c.dynamic_viscosity.value();
    let cross = c.bulk_viscosity.value() - (2.0 / 3.0) * mu;
    out.fill(0.0);
    out[n] = -longitudinal * normal_gradient;
    for gradient in transverse_gradients {
        out[n] -= cross * gradient;
    }
    let viscosity_scale = viscosity_scale::<TRANSVERSE>(longitudinal, cross, mu);
    let diffusivity_m2_s = viscosity_scale / p.density;
    if !out[n].is_finite() || !diffusivity_m2_s.is_finite() {
        return Err(DissipationError::Invalid("wall traction overflow"));
    }
    poll(cancelled)?;
    Ok(TransportFaceReport {
        diffusivity_m2_s,
        ..TransportFaceReport::default()
    })
}

// Conservative row-sum estimate for normal/transverse Newtonian couplings.
// A spatial owner supplies its own inverse lengths, summed over all axes.
fn viscosity_scale<const TRANSVERSE: usize>(longitudinal: f64, cross: f64, mu: f64) -> f64 {
    if TRANSVERSE == 0 {
        longitudinal
    } else {
        longitudinal + TRANSVERSE as f64 * cross.abs() + (2 * TRANSVERSE) as f64 * mu
    }
}
