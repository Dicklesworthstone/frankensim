//! Locally reacting compliant walls in the characteristic tube network.
//!
//! Reuse fs-phs's WallPin law: sigma x'' + r x' + K x = p, U = A x'.
//! Hence R = r/A, L = sigma/A, C = A/K. These are interior SHUNTS, not
//! terminal losses or output filters. Pressure drives wall motion, which returns
//! volume flow to the same junction solve. Mechanical storage and dissipated
//! wall work are included in the existing network's balance.
//!
//! This is a fixed-reference, linear locally reacting wall, NOT a shell model:
//! patches have no axial structural coupling or exterior radiation, and the
//! undeformed tube geometry remains fixed. The wall parameters must come from
//! the caller; a material name cannot select them. Small-displacement and spatial
//! refinement validity are the caller's responsibility. No gas viscothermal
//! boundary-layer law or automatic broadband identification is implied.

use super::network::{ApertureNetwork, NetworkNode, TubeNetworkSpec, TubeSection};
use super::tube::UniformTubeSpec;
use crate::acoustic_realize::AcousticRealizeError;
pub use fs_phs::WallPin;
use fs_vfit::impedance::SeriesImpedanceSpec;
use fs_vfit::relaxation::RelaxationImpedanceSpec;

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

/// Explicit area participating in the existing per-area wall law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallPatch {
    /// Wetted area [m^2], not the bore cross-sectional area.
    pub area_m2: f64,
    /// Surface mass [kg/m^2], stiffness [Pa/m], and resistance [Pa s/m].
    pub wall: WallPin,
}

/// Actual end-of-step wall motion, reconstructed from the accepted load state.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WallState {
    /// Displacement away from the bore [m]; positive increases its volume.
    pub displacement_m: f64,
    /// Outward velocity [m/s].
    pub velocity_m_s: f64,
    /// Independently reconstructed mechanical energy A(sigma v^2 + K x^2)/2 [J].
    pub stored_energy_j: f64,
}

impl WallPatch {
    /// Map the same fs-phs wall law into acoustic pressure/volume-flow units.
    ///
    /// # Errors
    /// Nonpositive mass, stiffness or area, negative resistance, or unrepresentable coefficients.
    pub fn impedance(&self) -> Result<SeriesImpedanceSpec, AcousticRealizeError> {
        let w = self.wall;
        if ![self.area_m2, w.surface_density, w.stiffness_per_area]
            .iter().all(|x| x.is_finite() && *x > 0.0)
            || !w.resistance.is_finite() || w.resistance < 0.0
        {
            return Err(invalid("wall patch requires positive finite area, surface mass and stiffness, and nonnegative resistance"));
        }
        let load = SeriesImpedanceSpec {
            resistance_pa_s_m3: w.resistance / self.area_m2,
            inertance_pa_s2_m3: w.surface_density / self.area_m2,
            compliance_m3_pa: Some(self.area_m2 / w.stiffness_per_area),
        };
        load.validate().map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        if load.inertance_pa_s2_m3 <= 0.0
            || (w.resistance > 0.0 && load.resistance_pa_s_m3 == 0.0)
        {
            return Err(invalid("wall patch conversion lost an explicit physical coefficient"));
        }
        Ok(load)
    }

    /// The interior shunt consumed directly by TubeNetworkSpec::nodes.
    ///
    /// # Errors
    /// Same physical admission as `impedance`.
    pub fn shunt(&self) -> Result<NetworkNode, AcousticRealizeError> {
        let load = RelaxationImpedanceSpec::new(self.impedance()?, &[])
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        Ok(NetworkNode::Shunt { load })
    }

    /// Read motion through this declared area/law map, requiring its acoustic
    /// impedance to match the node. This checks consistency, not unique material
    /// identification: different scaled geometries can have the same impedance.
    /// Never invent motion from an output-pressure waveform.
    ///
    /// # Errors
    /// Node/coefficient mismatch or nonfinite reconstructed coordinates/energy.
    pub fn observe(&self, network: &ApertureNetwork, node: usize) -> Result<WallState, AcousticRealizeError> {
        if network.spec().nodes.get(node).copied() != Some(self.shunt()?) {
            return Err(invalid("wall observation requires the same admitted patch at the named node"));
        }
        let state = network.load_state(node).ok_or_else(|| invalid("wall load state is unavailable"))?;
        let x = state.compliance_pressure_pa / self.wall.stiffness_per_area;
        let v = state.inertive_flow_m3_s / self.area_m2;
        let energy = 0.5 * self.area_m2 * (self.wall.surface_density * v * v
            + self.wall.stiffness_per_area * x * x);
        if ![x, v, energy].iter().all(|x| x.is_finite()) {
            return Err(invalid("wall motion reconstruction left the finite set"));
        }
        Ok(WallState { displacement_m: x, velocity_m_s: v, stored_energy_j: energy })
    }
}

/// A uniform tube with wall-area quadrature at the centers of equal axial strips.
/// Patch i belongs to node i+1. End half-sections preserve the requested total
/// length; patch areas cover the lateral wall once, without end discs or overlaps.
#[derive(Debug, Clone, PartialEq)]
pub struct LinedTube {
    /// Ready for the existing ApertureNetwork constructor.
    pub network: TubeNetworkSpec,
    /// Requested physical patch areas/law, in axial order.
    pub patches: Vec<WallPatch>,
}

/// Lower a cylindrical wall into independent, centered shunts and existing tube
/// sections. This is spatial lumping, not a new time integrator. Each half/end
/// section must still satisfy the runtime's positive integer-transit admission.
/// The total requested length-error allowance is divided among all sections;
/// patch areas use REQUESTED geometry and remain inspectable alongside the
/// runtime's represented section lengths. Refine space/time to assess this error.
///
/// `max_geometry_bytes` bounds requested builder vector payloads only. Runtime
/// propagation/load memory has the separate `tube.max_wave_memory_bytes` ceiling.
///
/// # Errors
/// Physical admission, zero/overflowing patch count, builder budget/allocation,
/// or unrepresentable strip lengths, areas or tolerances. No runtime is advanced.
pub fn lined_tube(
    tube: UniformTubeSpec, wall: WallPin, patch_count: usize, max_geometry_bytes: usize,
) -> Result<LinedTube, AcousticRealizeError> {
    if patch_count == 0 || ![tube.length_m, tube.radius_m, tube.sound_speed_m_s]
        .iter().all(|x| x.is_finite() && *x > 0.0)
        || !tube.terminal_reflection.is_finite() || tube.terminal_reflection.abs() > 1.0
        || !tube.max_length_error_m.is_finite() || tube.max_length_error_m < 0.0
    {
        return Err(invalid("lined tube requires positive geometry and patch count, a passive termination and finite length tolerance"));
    }
    let overflow = || invalid("lined tube geometry size overflow");
    let node_count = patch_count.checked_add(2).ok_or_else(overflow)?;
    let section_count = patch_count.checked_add(1).ok_or_else(overflow)?;
    let bytes = node_count.checked_mul(core::mem::size_of::<NetworkNode>())
        .and_then(|n| section_count.checked_mul(core::mem::size_of::<TubeSection>()).and_then(|s| n.checked_add(s)))
        .and_then(|n| patch_count.checked_mul(core::mem::size_of::<WallPatch>()).and_then(|s| n.checked_add(s)))
        .ok_or_else(overflow)?;
    if bytes > max_geometry_bytes { return Err(invalid("lined tube geometry payload exceeds budget")); }
    let dx = tube.length_m / patch_count as f64;
    let patch = WallPatch { area_m2: 2.0 * core::f64::consts::PI * tube.radius_m * dx, wall };
    let shunt = patch.shunt()?;
    let tolerance = tube.max_length_error_m / section_count as f64;
    if !dx.is_finite() || dx <= 0.0 || 0.5 * dx == 0.0
        || (tube.max_length_error_m > 0.0 && tolerance == 0.0)
    {
        return Err(invalid("lined tube strip length or tolerance is not representable"));
    }
    let mut nodes = Vec::new();
    let mut sections = Vec::new();
    let mut patches = Vec::new();
    nodes.try_reserve_exact(node_count).map_err(|_| invalid("lined tube allocation failed"))?;
    sections.try_reserve_exact(section_count).map_err(|_| invalid("lined tube allocation failed"))?;
    patches.try_reserve_exact(patch_count).map_err(|_| invalid("lined tube allocation failed"))?;
    nodes.push(NetworkNode::Inlet);
    for _ in 0..patch_count { nodes.push(shunt); patches.push(patch); }
    nodes.push(NetworkNode::Termination { reflection: tube.terminal_reflection });
    for i in 0..section_count {
        sections.push(TubeSection { nodes: [i, i + 1],
            length_m: if i == 0 || i == patch_count { 0.5 * dx } else { dx },
            radius_m: tube.radius_m, max_length_error_m: tolerance });
    }
    Ok(LinedTube { network: TubeNetworkSpec { nodes, sections,
        sound_speed_m_s: tube.sound_speed_m_s, max_wave_memory_bytes: tube.max_wave_memory_bytes }, patches })
}
