//! One material evaluation supplies electrical and fluid transport together.
//!
//! The existing current owner generates conductor-frame heat and Lorentz force;
//! the existing local-transport owner evolves momentum and sensible heat. The
//! initial-state coefficients and sources are held over one step. Endpoint
//! sampling checks material support, not a coupled midpoint or an induction law.

use super::{
    ChtError, ElectricHeatError, ElectricHeatInput, FluidTransportState,
    MagneticHeatStep, MaterialSlab, MovingSlab, ThermalLbm, ThermalLbm3,
    ThermalSourceError, VoxelDomain, admit, admit_field, balance, combined_forces,
    invalid, poll, positive, source_at, temperatures, velocities,
};
use crate::material::transport::{
    MaterialTransportError, TransportSample, candidate_with_samples, sample,
};
use fs_exec::CancelGate;
use fs_material::{MaterialError, fluid_transport::FluidTransport};

/// Electrical response and the concrete fluid-material owner at ONE query point.
/// Missing electrical data must return an error; `None` explicitly means an
/// insulator. Electrode cells still need fluid transport, but their resistivity
/// is ignored by the ideal-electrode model. Raw inputs do not acquire source
/// authority: retain the original source-backed object in `material`.
#[derive(Debug, Clone, PartialEq)]
pub struct ElectricalTransport<F> {
    material: F,
    resistivity_ohm_m: Option<f64>,
}
impl<F> ElectricalTransport<F> {
    /// Pair an existing material evaluation with isotropic Ohmic resistivity.
    /// Does not evaluate `material.transport_state()`; the shared sampler does
    /// so exactly once and retains that exact snapshot for all conversions.
    pub fn try_new(material: F, resistivity_ohm_m: Option<f64>) -> Result<Self, MaterialError> {
        if resistivity_ohm_m.is_some_and(|r| !r.is_finite() || r <= 0.0) {
            return Err(MaterialError::Parameters {
                what: "electrical transport needs positive finite resistivity or explicit insulation".into(),
            });
        }
        Ok(Self { material, resistivity_ohm_m })
    }
    /// Original concrete material object, including any retained source receipts.
    #[must_use]
    pub const fn material(&self) -> &F { &self.material }
    /// Isotropic resistivity [ohm m], or an explicit insulating declaration.
    #[must_use]
    pub const fn resistivity_ohm_m(&self) -> Option<f64> { self.resistivity_ohm_m }
}
impl<F: FluidTransport> FluidTransport for ElectricalTransport<F> {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.material.transport_state()
    }
}

/// Actual electromagnetic fields, thermal balance and both material epochs.
/// Electrical resistivity and local k/mu come from the SAME applied samples.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialMagneticStep<F> {
    /// Existing electrical/force/heat result. Mechanical work is still lagged
    /// source power times dt, not measured kinetic-energy change.
    pub step: MagneticHeatStep,
    /// Material objects and exact transport snapshots actually applied.
    /// Includes fluid electrode cells; excludes the thermal wall reservoirs.
    pub applied: Vec<TransportSample<ElectricalTransport<F>>>,
    /// Material evaluations at the accepted endpoint, not applied coefficients.
    pub endpoint: Vec<TransportSample<ElectricalTransport<F>>>,
}

/// No failure publishes a partial slab. Captured external callback state and
/// receipts, including those from rejected trials, remain caller-owned.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialMagneticError {
    /// Original material sampling/storage/transport failure, with cell context.
    Material(MaterialTransportError),
    /// Existing geometry, electrical, heat, transfer, velocity or cancellation error.
    Coupling(ElectricHeatError),
}
impl core::fmt::Display for MaterialMagneticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => write!(f, "{e}"),
            Self::Coupling(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for MaterialMagneticError {}
impl From<MaterialTransportError> for MaterialMagneticError {
    fn from(e: MaterialTransportError) -> Self { Self::Material(e) }
}
impl From<ElectricHeatError> for MaterialMagneticError {
    fn from(e: ElectricHeatError) -> Self { Self::Coupling(e) }
}
impl From<ChtError> for MaterialMagneticError {
    fn from(e: ChtError) -> Self { Self::Coupling(e.into()) }
}
impl From<ThermalSourceError> for MaterialMagneticError {
    fn from(e: ThermalSourceError) -> Self { Self::Coupling(e.into()) }
}

#[allow(clippy::too_many_arguments)]
fn advance<S: MovingSlab<D, Q>, F: FluidTransport, const D: usize, const Q: usize>(
    slab: &mut MaterialSlab<S>, domain: VoxelDomain, plane: usize,
    input: ElectricHeatInput<'_>, field_t: &[[f64; 3]], external: &[[f64; D]],
    depth_m: f64, isolated_current_tolerance_a: Option<f64>, gate: &CancelGate,
    mut resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
) -> Result<MaterialMagneticStep<F>, MaterialMagneticError> {
    poll(gate)?;
    let cells = domain.cell_count();
    admit_field::<D>(field_t, cells)?;
    if let Some(tol) = isolated_current_tolerance_a
        && (!tol.is_finite() || tol < 0.0 || !input.setup.electrodes.is_empty())
    {
        return Err(invalid("isolated material transport needs no electrodes and a finite nonnegative current tolerance").into());
    }
    let temperature = temperatures::<S, D, Q>(slab, input.limits.thermal, gate)?;
    let velocity = velocities::<S, D, Q>(slab, cells, plane, gate)?;
    let applied = sample::<S, ElectricalTransport<F>, D, Q>(
        slab, &slab.lattice, cells, plane, input.limits.thermal.temperature_range_k,
        gate, &mut resolve,
    )?;
    // The source solver requests only non-electrode interior cells. Do not
    // re-query a mutable material owner or use another cell/state's resistivity.
    let mut electrical = |cell: usize, _: f64, _: &FluidTransportState| {
        Ok(applied[cell - plane].material.resistivity_ohm_m())
    };
    let source = source_at(
        &domain, plane, input, field_t, &temperature, &velocity,
        &slab.force_units.transport, slab.force_units.timestep_s,
        isolated_current_tolerance_a, gate, &mut electrical,
    )?;
    let forces = combined_forces(external, &source.motion, gate)?;
    let (trial, heat) = candidate_with_samples::<S, ElectricalTransport<F>, D, Q>(
        slab, &source.power, &forces, depth_m, input.limits.thermal, &applied, gate,
    )?;
    let next = MaterialSlab {
        lattice: trial, top_temperature_k: slab.top_temperature_k,
        temperature_difference_k: slab.temperature_difference_k,
        force_units: slab.force_units,
    };
    velocities::<S, D, Q>(&next, cells, plane, gate)?;
    let endpoint = sample::<S, ElectricalTransport<F>, D, Q>(
        slab, &next.lattice, cells, plane, input.limits.thermal.temperature_range_k,
        gate, &mut resolve,
    )?;
    let step = balance(source, heat, input)?;
    poll(gate)?;
    *slab = next;
    Ok(MaterialMagneticStep { step, applied, endpoint })
}

impl MaterialSlab<ThermalLbm> {
    /// Advance wired 2D current, magnetic force, heat and LOCAL material transport.
    /// One callback at each accepted Kelvin point supplies both k/mu and Ohmic
    /// resistivity. Pressure, composition and source policy remain caller-owned.
    /// The callback includes fluid electrode cells (their electrical entry is
    /// unused), and is called again at each candidate endpoint before publish.
    ///
    /// B is in tesla and must be out-of-plane; zero B gives Ohmic heating with
    /// local transport. Forces are N/m3, depth is metres. Density/Cp must match
    /// retained storage; no mass or enthalpy is fabricated by coefficient edits.
    /// Rates, current and force are explicit initial-state values for one step.
    /// No induction, phase change, unconditional stability, total-energy or
    /// within-step coefficient convergence is claimed. Existing electrical
    /// boundaries, ideal electrodes and low-Mach conditions remain in force.
    /// All material, field, transport, endpoint and balance refusals are atomic.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_with_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 2]], depth_m: f64, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<F>, MaterialMagneticError> {
        self.step_material_motion(input, field_t, forces_n_m3, depth_m, None, gate, resolve)
    }

    /// Electrode-free counterpart, with an explicit local current tolerance [A].
    /// Reuses the isolated current owner; numerical gauges retain all material,
    /// heat and force contributions. All other semantics match the wired method.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_isolated_with_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 2]], depth_m: f64, current_tolerance_a: f64, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<F>, MaterialMagneticError> {
        self.step_material_motion(input, field_t, forces_n_m3, depth_m, Some(current_tolerance_a), gate, resolve)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_material_motion<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces: &[[f64; 2]], depth_m: f64, isolated: Option<f64>, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<F>, MaterialMagneticError> {
        let (_, plane) = admit::<ThermalLbm, 2, 5>(self, input, forces, 1, self.lattice.grid.ny, gate)?;
        let dx = self.force_units.cell_size_m;
        if !positive(depth_m) || !positive((dx * dx) * depth_m) {
            return Err(invalid("invalid represented 2D cell volume").into());
        }
        let domain = VoxelDomain::graded([
            vec![dx; self.lattice.grid.nx], vec![dx; self.lattice.grid.ny], vec![depth_m],
        ])?;
        advance::<ThermalLbm, F, 2, 5>(self, domain, plane, input, field_t, forces, depth_m, isolated, gate, resolve)
    }
}
impl MaterialSlab<ThermalLbm3> {
    /// 3D counterpart: retained dx^3 volumes, all vector components, D3Q7 heat
    /// scaling and the same single-evaluation material/circuit/transport path.
    pub fn step_magnetothermal_with_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 3]], gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<F>, MaterialMagneticError> {
        self.step_material_motion(input, field_t, forces_n_m3, None, gate, resolve)
    }

    /// 3D electrode-free counterpart, preserving the existing gauge/current
    /// continuity and conductor-frame heat semantics without phantom electrodes.
    #[allow(clippy::too_many_arguments)]
    pub fn step_magnetothermal_isolated_with_material_transport<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces_n_m3: &[[f64; 3]], current_tolerance_a: f64, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<F>, MaterialMagneticError> {
        self.step_material_motion(input, field_t, forces_n_m3, Some(current_tolerance_a), gate, resolve)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_material_motion<F: FluidTransport>(
        &mut self, input: ElectricHeatInput<'_>, field_t: &[[f64; 3]],
        forces: &[[f64; 3]], isolated: Option<f64>, gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<ElectricalTransport<F>, MaterialError>,
    ) -> Result<MaterialMagneticStep<F>, MaterialMagneticError> {
        let (_, plane) = admit::<ThermalLbm3, 3, 7>(self, input, forces, 2, self.lattice.grid.nz, gate)?;
        let g = &self.lattice.grid;
        let dx = self.force_units.cell_size_m;
        if !positive((dx * dx) * dx) { return Err(invalid("invalid represented 3D cell volume").into()); }
        let domain = VoxelDomain::new(g.nx, g.ny, g.nz, dx)?;
        advance::<ThermalLbm3, F, 3, 7>(self, domain, plane, input, field_t, forces, dx, isolated, gate, resolve)
    }
}
