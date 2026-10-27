//! Distributed interior air on the SAME material-shell modes as every string.
//!
//! The original tetrahedral Helmholtz owner supplies the pressure basis. Its
//! boundary quadrature integrates each original P1 shell normal trace against
//! each pressure shape. The existing CavityCoupling then supplies reciprocal
//! compression/standing-wave storage to the existing complete friction solve.
//! No assigned body frequencies, fitted timbre, independent air stepper, or
//! duplicated uniform compliance is introduced.
//!
//! Scope: one connected conforming volume with an explicit partition between
//! material shell faces and optional fixed compact fluid mouths, frozen coherent
//! gas, linear shell kinematics and optional source-gas heat diffusion. Mesh, pressure
//! order, retained frequency window and structural modes require convergence.
//! A declared mouth supplies its internal impedance; the existing exterior BEM
//! supplies its joint pressure reaction. This does not resolve a distributed
//! f-hole field, supply viscous wall shear, or claim instrument validation.
//! Thermal loss requires explicitly isothermal solid walls and a separately
//! admitted heat FE/eigenwindow; the default remains lossless interior air.

use super::*;
use crate::render::plate::impact::{ImpactError, cavity::{
    CavityCoupling, MAX_CAVITY_PRESSURE_MODES,
}};
use crate::shear_shell::ResolvedShearShell;
use crate::vibroacoustic::CavityModes;
use fs_blake3::ContentHash;
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, HelmholtzModes};
use fs_math::det;
use fs_modal::SliceOptions;
use fs_plate::shell::pressure::ShellFace;
use std::collections::BTreeMap;

mod openings;
pub use openings::{MAX_SHELL_CAVITY_OPENINGS, ShellCavityOpening};
mod thermal;
pub use thermal::{ShellCavityThermalLoss, ShellCavityThermalOptions};

/// Numerical admission, independent of instrument name or material choice.
#[derive(Clone, Copy, Debug)]
pub struct ShellCavityOptions {
    /// Positive retained pressure-frequency ceiling [Hz]. Uniform compression
    /// is always included; every mode inside the window must fit the allowance.
    pub maximum_hz: f64,
    pub eigensolver: SliceOptions,
    /// 1..=256 complete pressure modes, including the uniform mode.
    pub maximum_pressure_modes: usize,
    /// 4..=8192 original pressure degrees of freedom.
    pub maximum_dofs: usize,
    /// Complete selected shell/volume boundary, 4..=8192 triangles.
    pub maximum_boundary_faces: usize,
    /// Face-quadrature and complete modal projection work, bounded before solving.
    pub maximum_projection_terms: usize,
    /// Original per-component acoustic time owner admission.
    pub acoustic_budget: ModalAcousticTimeBudget,
}

/// Actual source geometry, pressure basis and reciprocal storage retained by
/// single/multiple bows, material-memory rebuilds and scheduled performance.
#[derive(Clone, Debug)]
pub struct MaterialShellCavity {
    gas: GasState,
    basis_identity: ContentHash,
    positions: Vec<[f64; 3]>,
    tetrahedra: Vec<[u32; 4]>,
    order: usize,
    options: ShellCavityOptions,
    modal_report: HelmholtzModes,
    pressure_shapes: Vec<Vec<f64>>,
    overlap_rows: Vec<f64>,
    volume_m3: f64,
    coupling: CavityCoupling,
    openings: Vec<openings::PreparedOpening>,
    thermal_loss: Option<ShellCavityThermalLoss>,
}
impl MaterialShellCavity {
    pub const fn gas(&self) -> &GasState { &self.gas }
    pub const fn basis_identity(&self) -> ContentHash { self.basis_identity }
    pub fn positions(&self) -> &[[f64; 3]] { &self.positions }
    pub fn tetrahedra(&self) -> &[[u32; 4]] { &self.tetrahedra }
    pub const fn order(&self) -> usize { self.order }
    pub const fn options(&self) -> &ShellCavityOptions { &self.options }
    /// Original eigensolver report, before replacing its roundoff-level uniform
    /// eigenvector by the exact normalized constant of the original FE space.
    pub const fn modal_report(&self) -> &HelmholtzModes { &self.modal_report }
    /// Actual pressure shapes in original FE degree-of-freedom order. Shape 0
    /// is the exact normalized constant; all other original modes are retained.
    pub fn pressure_shapes(&self) -> &[Vec<f64>] { &self.pressure_shapes }
    /// Row-major structure-by-pressure work integral, with outward fluid signs.
    pub fn overlap_rows(&self) -> &[f64] { &self.overlap_rows }
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    pub const fn coupling(&self) -> &CavityCoupling { &self.coupling }
    /// Original gas-conduction reduction, when explicitly selected. With this
    /// law, use the bank's pressure observer, which includes accepted thermal
    /// memory; `coupling()` alone describes only the adiabatic storage basis.
    pub const fn thermal_loss(&self) -> Option<&ShellCavityThermalLoss> {
        self.thermal_loss.as_ref()
    }
    /// Caller-owned interleaved displacement/velocity scratch for observation.
    pub fn observation_state_len(&self) -> usize { 2 * self.coupling.total_modes() }
    /// Actual compact fluid apertures, in the source's declared order.
    pub fn opening_count(&self) -> usize { self.openings.len() }
    pub fn opening(&self, index: usize) -> Option<&ShellCavityOpening> {
        self.openings.get(index).map(|opening| &opening.source)
    }
    /// Bounding mouth radius about its area centroid or internal length,
    /// whichever is larger; the compact chart bounds the whole physical patch.
    pub fn opening_extent_m(&self, index: usize) -> Option<f64> {
        self.openings.get(index).map(|opening| opening.extent_m)
    }
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn opening_triangles(
        &self, index: usize,
    ) -> Option<&[[[f64; 3]; 3]]> {
        self.openings.get(index).map(|opening| opening.triangles.as_slice())
    }
    /// Body modes followed by true neck coordinates, relative to the original
    /// body component. Interior standing-wave inertias are never exterior ports.
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn exterior_body_coordinates(
        &self,
    ) -> Result<Vec<usize>, ShellCavityError> {
        let mut coordinates: Vec<_> = (0..self.coupling.structural_modes()).collect();
        for index in 0..self.openings.len() {
            coordinates.push(self.coupling.neck_radiation_port(index)
                .map_err(ShellCavityError::Coupling)?.coordinate);
        }
        self.coupling.admit_exterior_sources(&coordinates).map_err(ShellCavityError::Coupling)?;
        Ok(coordinates)
    }
}

#[derive(Debug)]
pub enum ShellCavityError {
    Invalid(&'static str),
    Acoustic(AcousticFeError),
    Coupling(ImpactError),
    Bank(ReciprocalStringError),
    Radiation(super::radiation::bem::ShellBemError),
    Cancelled,
}
impl core::fmt::Display for ShellCavityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material shell cavity: {self:?}")
    }
}
impl std::error::Error for ShellCavityError {}
fn check_cancel(gate: &CancelGate) -> Result<(), ShellCavityError> {
    if gate.is_requested() { Err(ShellCavityError::Cancelled) } else { Ok(()) }
}
fn bad(what: &'static str) -> ShellCavityError { ShellCavityError::Invalid(what) }

impl MaterialStringBankSource {
    pub const fn spatial_cavity(&self) -> Option<&MaterialShellCavity> {
        self.spatial_cavity.as_ref()
    }
}

impl ReciprocalStringBank {
    /// Join the original 3D pressure FE to the actual material shell or assembled
    /// material body. All strings, bridge/torque reactions, full body viscosity,
    /// body observation, original material histories and exterior loading remain
    /// in the same reciprocal owner. A later string gas-loading or torque rebuild
    /// retains the added air component and every compression connection.
    ///
    /// Every cavity boundary face must coincide EXACTLY with one explicitly
    /// selected original shell cavity facet, with the same outward-from-gas
    /// orientation. No welding tolerance, hole filling or nearest-face guessing.
    /// The FE medium must be the bank's coherent gas density and sound speed.
    /// Its volume, order and pressure window are data; material substitution
    /// rebuilds the original shell modes and projects the same volume anew.
    ///
    /// Cold construction only, before string/fluid motion or nonlinear binding.
    /// Existing body initial displacement/velocity and material memory survive;
    /// an initial displaced wall really compresses the air. Interior pressure
    /// changes exterior sound through the body's accepted loaded acceleration.
    /// Interior pressure coordinates are not relabeled as microphone channels.
    pub fn with_material_shell_cavity(
        self, fe: &HelmholtzFe<'_>, options: ShellCavityOptions, gate: &CancelGate,
    ) -> Result<Self, ShellCavityError> {
        self.with_material_shell_cavity_inner(fe, options, Vec::new(), None, gate)
    }

    /// Connect explicitly declared compact mouths to the same original FE and
    /// body. Every mouth replaces actual missing solid faces, not a rigid cap.
    /// Its area and pressure averages come from the original FE face integrals;
    /// only internal neck length and resistance are supplied. Mouth rims must
    /// be fixed in the retained structural basis. All initial neck states are
    /// zero; normal network preload can subsequently excite them.
    ///
    /// This returns the explicit zero-exterior-pressure comparison. Attach
    /// `with_material_shell_bem` for one joint body/mouth load, or use the atomic
    /// constructor below. No exterior inertia or loss is counted in the neck.
    pub fn with_material_shell_cavity_openings(
        self, fe: &HelmholtzFe<'_>, options: ShellCavityOptions,
        openings: Vec<ShellCavityOpening>, gate: &CancelGate,
    ) -> Result<Self, ShellCavityError> {
        self.with_material_shell_cavity_inner(fe, options, openings, None, gate)
    }

    /// Couple gas heat diffusion on the SAME original cavity FE to every solid,
    /// standing-wave and compact-mouth compression coordinate. Solid boundary
    /// faces are explicitly isothermal; apertures have zero diffusive heat flux.
    /// All thermal states join the original mechanical/material propagator.
    /// Select optional body/string constitutive memories before this cold call.
    /// Their histories remain an unchanged prefix; new gas history starts
    /// unrelaxed, preserving the original instantaneous adiabatic pressure.
    /// This is thermal conduction only, not viscous wall shear or finite-wall
    /// heat storage. The declared missing-projection bound is binding.
    pub fn with_material_shell_cavity_thermal(
        self, fe: &HelmholtzFe<'_>, options: ShellCavityOptions,
        openings: Vec<ShellCavityOpening>, thermal: ShellCavityThermalOptions,
        gate: &CancelGate,
    ) -> Result<Self, ShellCavityError> {
        thermal.validate()?;
        self.with_material_shell_cavity_inner(fe, options, openings, Some(thermal), gate)
    }

    /// One source-derived interior/compact-mouth/exterior construction. The
    /// existing BEM and passive realization retain all signed body-mouth and
    /// mouth-mouth reactions in the original exact mechanical propagator.
    /// Independent columns use the original symmetric P1/P0 source and must
    /// pass unchanged reciprocity, passivity and causal-fit gates. Its geometric
    /// resonance-exclusion and angular-resolution chart remains authoritative.
    pub fn with_material_shell_cavity_and_bem(
        self, fe: &HelmholtzFe<'_>, options: ShellCavityOptions,
        openings: Vec<ShellCavityOpening>,
        exterior: super::radiation::bem::ShellBemOptions, gate: &CancelGate,
    ) -> Result<(Self, super::radiation::bem::ShellBemReport), ShellCavityError> {
        self.with_material_shell_cavity_inner(fe, options, openings, None, gate)?
            .with_material_shell_bem(exterior, gate).map_err(ShellCavityError::Radiation)
    }

    fn with_material_shell_cavity_inner(
        mut self, fe: &HelmholtzFe<'_>, options: ShellCavityOptions,
        opening_sources: Vec<ShellCavityOpening>, thermal: Option<ShellCavityThermalOptions>,
        gate: &CancelGate,
    ) -> Result<Self, ShellCavityError> {
        check_cancel(gate)?;
        let string_count = self.source.strings.len();
        let body_count = self.source.fluid_loading.as_ref().map_or(
            self.system.components().len() - string_count, |f| f.body_components,
        );
        if self.source.spatial_cavity.is_some() || self.source.material_cavity.is_some()
            || (!opening_sources.is_empty() && self.source.body_radiation.is_some())
            || self.system.samples_rendered() != 0 || self.system.has_nonlinear_potentials()
            || body_count != 1
            || self.system.components().iter().enumerate()
                .filter(|(i,_)| *i != string_count)
                .flat_map(|(_,c)| c.states()).any(|s|
                    s.displacement_m_sqrt_kg != 0.0 || s.velocity_m_sqrt_kg_per_s != 0.0)
            || self.system.radiation_state().iter().any(|v| *v != 0.0)
        {
            return Err(bad("distributed cavity requires one cold material body without duplicated air or discarded motion"));
        }
        self.source.gas.validate_coherence().map_err(|_| bad("incoherent body cavity gas"))?;
        if fe.medium().density_kg_m3.to_bits() != self.source.gas.density.to_bits()
            || fe.medium().sound_speed_m_s.to_bits() != self.source.gas.sound_speed.to_bits()
        {
            return Err(bad("cavity pressure basis and original material bank must share the same coherent gas"));
        }
        let body = self.system.components()[string_count].clone();
        let (identity, _) = self.source.body_radiation_basis(&self.system, gate)
            .map_err(ShellCavityError::Bank)?;
        let (binding, area) = prepare(&self.source, fe, options, identity,
            body.modes().len(), self.source.options.network.max_modes, opening_sources, thermal, gate)?;
        let mut image = if binding.openings.is_empty() {
            binding.coupling.modal_network(body, self.source.options.sample_rate_hz,
                area, options.acoustic_budget, gate)
        } else {
            binding.coupling.modal_network_for_exterior_load(body,
                self.source.options.sample_rate_hz, area, options.acoustic_budget,
                &binding.exterior_body_coordinates()?, gate)
        }.map_err(ShellCavityError::Coupling)?;
        if let Some(thermal) = &binding.thermal_loss {
            thermal.replace_compression_links(&mut image, &binding.coupling, gate)?;
        }
        let added_components = image.components.len() - 1;
        let added_modes = binding.coupling.total_modes() - binding.coupling.structural_modes();
        let mut links = self.source.body_connections.clone();
        links.extend(image.connections);
        let mut fluid = self.source.fluid_loading.clone();
        if let Some(loading) = &mut fluid { loading.body_components += added_components; }
        let bridge_rows: Vec<_> = self.source.strings.iter().map(|s| s.connections.as_slice()).collect();
        let mut system = assemble_sources(&self.cards, Some(&self.source.strings),
            image.components, &bridge_rows, &links, &self.source.translational_body_connections,
            &self.source.rotational_body_connections,
            fluid.as_ref(), self.source.prepared_linear, &self.source.constant_loss_bridges(),
            self.source.residual_flexibility, self.system.port().law.clone(), &self.source.options, gate).map_err(ShellCavityError::Bank)?;
        if let Some(original) = self.system.acceleration_pressure_weights() {
            let insert = self.system.components()[..string_count+1].iter()
                .map(|component| component.modes().len()).sum::<usize>();
            let mut weights = Vec::with_capacity(system.mode_count());
            weights.extend_from_slice(&original[..insert]);
            weights.resize(insert + added_modes, 0.0);
            weights.extend_from_slice(&original[insert..]);
            system = system.with_acceleration_pressure(weights)
                .map_err(|e| ShellCavityError::Bank(ReciprocalStringError::Friction(e)))?;
        }
        if let Some(interface) = self.system.interface() {
            system = system.with_interface(interface.clone())
                .map_err(|e| ShellCavityError::Bank(ReciprocalStringError::Friction(e)))?;
        }
        // The source must describe the new complete physical prefix before
        // material-interface frequency checks or radiation restoration run.
        self.source.body_connections = links;
        self.source.fluid_loading = fluid;
        let mut memory = self.system.material_memory().to_vec();
        if let Some(thermal) = &binding.thermal_loss {
            memory.resize(memory.len()+thermal.memory_states(), 0.0);
        }
        self.source.spatial_cavity = Some(binding);
        system = self.source.restore_material_memory(system, &self.cards,
            Some(memory), gate).map_err(ShellCavityError::Bank)?;
        self.source.validate_body_memory_frequencies(&system, gate).map_err(ShellCavityError::Bank)?;
        self.source.validate_material_bridge_frequencies(&system, gate).map_err(ShellCavityError::Bank)?;
        self.system = system;
        check_cancel(gate)?;
        Ok(self)
    }

    /// Interior pressure coefficients [Pa] at the accepted complete bank state.
    /// Both buffers are caller-owned; no allocation or state advance occurs.
    pub fn spatial_cavity_pressures_into(&self, scratch: &mut [f64], pressure: &mut [f64])
        -> Result<(), ShellCavityError>
    {
        observe(&self.source, self.system.components(), self.system.material_memory(), scratch, pressure)
    }

    /// Accepted internal pressure, volume flow and kinetic storage of the
    /// original neck. This readout adds no force, energy or time integration.
    pub fn spatial_cavity_neck_observation(&self, index: usize, scratch: &mut [f64])
        -> Result<crate::render::plate::impact::cavity::neck::NeckObservation, ShellCavityError>
    {
        observe_neck(&self.source, self.system.components(), self.system.material_memory(), index, scratch)
    }

    /// Mouth-averaged exterior load pressure [Pa], positive against outward
    /// flow. It is already included in the original mechanical force and work
    /// balance. Zero denotes the explicitly unloaded comparison, not radiation.
    pub fn spatial_cavity_neck_exterior_pressure(&self, index: usize)
        -> Result<f64, ShellCavityError>
    {
        let cavity = self.source.spatial_cavity.as_ref()
            .ok_or_else(||bad("bank has no distributed cavity"))?;
        let port = cavity.coupling.neck_radiation_port(index).map_err(ShellCavityError::Coupling)?;
        let Some(load) = self.source.body_radiation.as_ref() else { return Ok(0.0); };
        let column = load.source_body_coordinates().iter().position(|&c|c == port.coordinate)
            .ok_or_else(||bad("exterior load lost its actual neck source"))?;
        let state = self.system.radiation_state();
        if state.is_empty() && load.model().poles.iter()
            .all(|pole| pole.coupling.iter().all(|&value| value == 0.0))
        { return Ok(0.0); }
        if state.len() != 2*load.model().poles.len() {
            return Err(bad("exterior neck pressure lost its original acoustic history"));
        }
        let pressure = load.model().poles.iter().zip(state.chunks_exact(2))
            .map(|(pole,state)|pole.coupling[column]*state[1]).sum::<f64>()
            /port.volume_weight_m2_per_sqrt_kg;
        if !pressure.is_finite() { return Err(bad("exterior neck pressure overflow")); }
        Ok(pressure)
    }
}
impl multiple::MultiBowStringBank {
    /// The same interior pressure after simultaneous bow/contact acceptance.
    pub fn spatial_cavity_pressures_into(&self, scratch: &mut [f64], pressure: &mut [f64])
        -> Result<(), ShellCavityError>
    {
        observe(self.source().bank(), self.system().components(), self.system().material_memory(), scratch, pressure)
    }

    /// The same accepted physical neck after simultaneous bow/contact solving.
    pub fn spatial_cavity_neck_observation(&self, index: usize, scratch: &mut [f64])
        -> Result<crate::render::plate::impact::cavity::neck::NeckObservation, ShellCavityError>
    {
        observe_neck(self.source().bank(), self.system().components(), self.system().material_memory(), index, scratch)
    }
}
fn observe(source: &MaterialStringBankSource, components: &[ModalAcousticTimeModel],
    memory: &[f64], scratch: &mut [f64], pressure: &mut [f64]) -> Result<(), ShellCavityError>
{
    let cavity = observation_state(source, components, scratch)?;
    if pressure.len() != cavity.coupling.cavity_modes() {
        return Err(bad("cavity observation requires the complete pressure buffer"));
    }
    // Stage the complete physical readout so an invalid memory suffix cannot
    // publish an adiabatic-only pressure as if it were the accepted endpoint.
    let mut values = [0.0; MAX_CAVITY_PRESSURE_MODES];
    let values = &mut values[..pressure.len()];
    cavity.coupling.adiabatic_pressures_into(scratch, values).map_err(ShellCavityError::Coupling)?;
    if let Some(thermal) = &cavity.thermal_loss { thermal.correct_pressure(memory, values)?; }
    pressure.copy_from_slice(values);
    Ok(())
}

fn observe_neck(source: &MaterialStringBankSource, components: &[ModalAcousticTimeModel],
    memory: &[f64], index: usize, scratch: &mut [f64])
    -> Result<crate::render::plate::impact::cavity::neck::NeckObservation, ShellCavityError>
{
    let cavity = observation_state(source, components, scratch)?;
    let mut pressure = [0.0; MAX_CAVITY_PRESSURE_MODES];
    let pressure = &mut pressure[..cavity.coupling.cavity_modes()];
    cavity.coupling.adiabatic_pressures_into(scratch,pressure).map_err(ShellCavityError::Coupling)?;
    if let Some(thermal) = &cavity.thermal_loss { thermal.correct_pressure(memory,pressure)?; }
    cavity.coupling.neck_observation_with_pressures(scratch,index,pressure).map_err(ShellCavityError::Coupling)
}
fn observation_state<'a>(source: &'a MaterialStringBankSource,
    components: &[ModalAcousticTimeModel], scratch: &mut [f64])
    -> Result<&'a MaterialShellCavity, ShellCavityError>
{
    let cavity = source.spatial_cavity.as_ref().ok_or_else(||bad("bank has no distributed cavity"))?;
    if scratch.len() != cavity.observation_state_len() {
        return Err(bad("cavity observation requires the complete declared state buffer"));
    }
    let start = source.strings.len();
    let count = 1 + usize::from(cavity.coupling.total_modes() > cavity.coupling.structural_modes());
    let bodies = components.get(start..start+count)
        .ok_or_else(|| bad("cavity observation lost its original body/air components"))?;
    if bodies[0].modes().len() != cavity.coupling.structural_modes()
        || bodies.iter().map(|c| c.states().len()).sum::<usize>() != cavity.coupling.total_modes()
    { return Err(bad("cavity observation lost its complete original modal basis")); }
    for (state, pair) in bodies.iter().flat_map(|body| body.states()).zip(scratch.chunks_exact_mut(2)) {
        pair[0] = state.displacement_m_sqrt_kg;
        pair[1] = state.velocity_m_sqrt_kg_per_s;
    }
    Ok(cavity)
}

type ShellPart<'a> = (&'a ResolvedShearShell, &'a [Vec<f64>]);
type Trace = ([u32; 3], [f64; 3]);
type FaceKey = [[u64; 3]; 3];
fn key(points: [[f64; 3]; 3]) -> FaceKey {
    let mut key = points.map(|p| p.map(|v| if v == 0.0 { 0 } else { v.to_bits() }));
    key.sort_unstable();
    key
}
fn parts(source: &MaterialStringBankSource) -> Result<Vec<ShellPart<'_>>, ShellCavityError> {
    if let Some(image) = &source.material_shell_image {
        Ok(vec![(image.source(), image.physical_modes())])
    } else if let Some(image) = &source.material_shell_assembly_image {
        Ok(image.source().components.iter().zip(image.physical_modes())
            .map(|(part, rows)| (&part.specimen, rows.as_slice())).collect())
    } else { Err(bad("distributed cavity requires the retained original material-shell basis")) }
}

fn prepare(source: &MaterialStringBankSource, fe: &HelmholtzFe<'_>, options: ShellCavityOptions,
    identity: ContentHash, body_modes: usize, maximum_mechanical_modes: usize,
    opening_sources: Vec<ShellCavityOpening>, thermal: Option<ShellCavityThermalOptions>, gate: &CancelGate)
    -> Result<(MaterialShellCavity, f64), ShellCavityError>
{
    let parts = parts(source)?;
    let boundary_count = fe.boundary_face_count();
    if !options.maximum_hz.is_finite() || options.maximum_hz <= 0.0
        || !(1..=MAX_CAVITY_PRESSURE_MODES).contains(&options.maximum_pressure_modes)
        || !(4..=8192).contains(&options.maximum_dofs) || fe.ndof() > options.maximum_dofs
        || !(4..=8192).contains(&options.maximum_boundary_faces)
        || boundary_count > options.maximum_boundary_faces || boundary_count < 4
        || options.maximum_projection_terms == 0 || options.maximum_projection_terms > 1_000_000_000
        || core::f64::consts::TAU * options.maximum_hz
            >= 0.9 * core::f64::consts::PI * f64::from(source.options.sample_rate_hz)
    { return Err(bad("complete cavity geometry, frequency window or pressure budget is not admitted")); }
    let order = fe.order();
    let local = (order+1)*(order+2)*(order+3)/6;
    let work = boundary_count.checked_mul(local).and_then(|n| n.checked_mul((order+1)*(order+1)))
        .and_then(|n| fe.ndof().checked_mul(options.maximum_pressure_modes).and_then(|p| n.checked_add(p)))
        .and_then(|n| n.checked_mul(body_modes.checked_add(opening_sources.len())?));
    let mass_work = fe.tetrahedra().len().checked_mul(local).and_then(|n| n.checked_mul(local))
        .and_then(|n| n.checked_mul(3 * options.maximum_pressure_modes + 4));
    let work = work.zip(mass_work).and_then(|(a,b)| a.checked_add(b));
    if work.is_none_or(|n| n > options.maximum_projection_terms) {
        return Err(bad("complete structural/acoustic overlap exceeds its projection allowance"));
    }
    connected_volume(fe, gate)?;
    let boundary = fe.boundary_faces();
    let openings = openings::prepare(&parts, fe, opening_sources, gate)?;
    let (traces, area) = surface_traces(&parts, fe, &boundary, body_modes, &openings, gate)?;
    let modal_report = fe.modes((-1.0, options.maximum_hz), &options.eigensolver)
        .map_err(ShellCavityError::Acoustic)?;
    check_cancel(gate)?;
    let count = modal_report.shapes.len();
    if count == 0 || count > options.maximum_pressure_modes
        || modal_report.below_low != 0 || modal_report.below_high != count
        || modal_report.frequency_hz.len() != count
        || modal_report.lambda.len() != count
    { return Err(bad("entire cavity eigenwindow, including uniform compression, must be retained")); }
    let (constant, volume) = constant_mode(fe, &modal_report.shapes[0], gate)?;
    let mut pressure_shapes = modal_report.shapes.clone();
    pressure_shapes[0] = constant;
    let mut frequencies = Vec::with_capacity(count);
    let mut norms = Vec::with_capacity(count);
    for (index, shape) in pressure_shapes.iter().enumerate() {
        check_cancel(gate)?;
        if shape.len() != fe.ndof() || shape.iter().any(|v| !v.is_finite())
            || (index != 0 && (!modal_report.lambda[index].is_finite() || modal_report.lambda[index] <= 0.0))
        { return Err(bad("cavity requires one constant mode and finite positive standing-wave eigenpairs")); }
        let norm = mass_product(fe, shape, shape, gate)?;
        if !norm.is_finite() || norm <= 0.0 {
            return Err(bad("cavity pressure norm is unrepresentable"));
        }
        // The exact uniform mode replaces only the known Neumann nullspace,
        // never a fitted frequency. Reject a disconnected/additional nullspace.
        if index != 0 && mass_product(fe, &pressure_shapes[0], shape, gate)?.abs()
            > 1e-7 * det::sqrt(norm)
        { return Err(bad("pressure eigenbasis does not separate uniform compression from standing waves")); }
        frequencies.push(if index == 0 { 0.0 }
            else { core::f64::consts::TAU * modal_report.frequency_hz[index] });
        norms.push(norm);
    }
    let mut overlap_rows = Vec::with_capacity(body_modes * count);
    for trace in &traces {
        check_cancel(gate)?;
        let load = fe.boundary_linear_load(trace).map_err(ShellCavityError::Acoustic)?;
        for shape in &pressure_shapes {
            let value = load.iter().zip(shape).map(|(l,p)| l*p).sum::<f64>();
            if !value.is_finite() { return Err(bad("cavity structural overlap overflow")); }
            overlap_rows.push(value);
        }
    }
    let basis = CavityModes { omegas: frequencies, lambdas: norms,
        interface: vec![Vec::new(); count], loss_factor: 0.0,
        rho0: source.gas.density, c0: source.gas.sound_speed };
    let mut coupling = CavityCoupling::new_with_mode_budget(&basis, body_modes,
        &overlap_rows, &vec![0.0; count], maximum_mechanical_modes)
        .map_err(ShellCavityError::Coupling)?;
    if !openings.is_empty() {
        let necks = openings.iter().map(|opening| opening.neck(fe, &pressure_shapes, gate))
            .collect::<Result<Vec<_>, _>>()?;
        coupling = coupling.with_necks_for_exterior_load(necks, gate)
            .map_err(ShellCavityError::Coupling)?;
        for (index,opening) in openings.iter().enumerate() {
            let omega = coupling.neck_fixed_wall_angular_frequency_rad_s(index)
                .map_err(ShellCavityError::Coupling)?
                .max(core::f64::consts::TAU*options.maximum_hz);
            if omega*opening.extent_m/source.gas.sound_speed > 0.3 {
                return Err(bad("complete mouth patch exceeds the compact acoustic chart"));
            }
        }
    }
    check_cancel(gate)?;
    let thermal_loss = thermal.map(|thermal| ShellCavityThermalLoss::prepare(fe,
        &source.gas, &pressure_shapes, &coupling, &openings, thermal, &options.eigensolver, gate))
        .transpose()?;
    if thermal_loss.is_some() { coupling.require_compression_memory(); }
    Ok((MaterialShellCavity { gas: source.gas, basis_identity: identity,
        positions: fe.positions().to_vec(), tetrahedra: fe.tetrahedra().to_vec(),
        order, options, modal_report, pressure_shapes, overlap_rows, volume_m3: volume,
        coupling, openings, thermal_loss }, area))
}

fn connected_volume(fe: &HelmholtzFe<'_>, gate: &CancelGate) -> Result<(), ShellCavityError> {
    fn root(parents: &mut [usize], mut node: usize) -> usize {
        while parents[node] != node {
            parents[node] = parents[parents[node]];
            node = parents[node];
        }
        node
    }
    let mut parents: Vec<_> = (0..fe.tetrahedra().len()).collect();
    let mut owners = BTreeMap::<[u32; 3], (usize, usize)>::new();
    for (cell, tet) in fe.tetrahedra().iter().enumerate() {
        check_cancel(gate)?;
        for omit in 0..4 {
            let mut face = [0; 3];
            let mut next = 0;
            for (index, &vertex) in tet.iter().enumerate() {
                if index != omit { face[next] = vertex; next += 1; }
            }
            face.sort_unstable();
            if let Some((other, count)) = owners.get_mut(&face) {
                *count += 1;
                if *count > 2 { return Err(bad("cavity volume has a nonmanifold internal face")); }
                let a = root(&mut parents, cell);
                let b = root(&mut parents, *other);
                parents[a] = b;
            } else { owners.insert(face, (cell, 1)); }
        }
    }
    if parents.is_empty() { return Err(bad("cavity has no volume cells")); }
    let first = root(&mut parents, 0);
    for cell in 1..parents.len() {
        if root(&mut parents, cell) != first {
            return Err(bad("one distributed cavity must be connected through complete volume faces"));
        }
    }
    Ok(())
}

fn mass_product(fe: &HelmholtzFe<'_>, a: &[f64], b: &[f64], gate: &CancelGate)
    -> Result<f64, ShellCavityError>
{
    if a.len() != fe.ndof() || b.len() != fe.ndof() {
        return Err(bad("pressure normalization lost its original FE dimensions"));
    }
    let mut sum = 0.0;
    for (i, value) in a.iter().enumerate() {
        if i % 256 == 0 { check_cancel(gate)?; }
        let (columns, entries) = fe.mass().row(i);
        sum += value * columns.iter().zip(entries).map(|(j,m)| m*b[*j]).sum::<f64>();
    }
    Ok(sum)
}
fn constant_mode(fe: &HelmholtzFe<'_>, eigenvector: &[f64], gate: &CancelGate)
    -> Result<(Vec<f64>, f64), ShellCavityError>
{
    // Hierarchical H1 has its vertex functions first. Their partition of unity
    // is the exact constant; every edge/face/cell-bubble coefficient is zero.
    let mut constant = vec![0.0; fe.ndof()];
    constant[..fe.positions().len()].fill(1.0);
    let volume = mass_product(fe, &constant, &constant, gate)?;
    if !volume.is_finite() || volume <= 0.0 { return Err(bad("invalid cavity volume")); }
    let scale = 1.0 / det::sqrt(volume);
    for value in &mut constant { *value *= scale; }
    let sign = if mass_product(fe, &constant, eigenvector, gate)? < 0.0 { -1.0 } else { 1.0 };
    let difference: Vec<_> = eigenvector.iter().zip(&constant).map(|(a,b)| a-sign*b).collect();
    let defect = mass_product(fe, &difference, &difference, gate)?;
    if !defect.is_finite() || defect > 1e-12 {
        return Err(bad("first pressure eigenvector is not the connected cavity's uniform mode"));
    }
    Ok((constant, volume))
}

fn surface_traces(parts: &[ShellPart<'_>], fe: &HelmholtzFe<'_>, boundary: &[[u32; 3]],
    modes: usize, openings: &[openings::PreparedOpening], gate: &CancelGate)
    -> Result<(Vec<Vec<Trace>>, f64), ShellCavityError>
{
    let mut faces = BTreeMap::new();
    let mut area = 0.0;
    for &face in boundary {
        check_cancel(gate)?;
        let (points, normal, face_area) = fe.boundary_geometry(face).map_err(ShellCavityError::Acoustic)?;
        if faces.insert(key(points), (face, points, normal)).is_some() {
            return Err(bad("cavity boundary has geometrically duplicate faces"));
        }
        area += face_area;
    }
    let mut traces = vec![Vec::with_capacity(boundary.len()); modes];
    for &(shell, rows) in parts {
        check_cancel(gate)?;
        let mesh = shell.mesh();
        let surfaces = shell.pressure_surfaces()
            .ok_or_else(|| bad("material shell must explicitly declare its interior pressure ownership"))?;
        if rows.len() != modes || rows.iter().any(|row|
            row.len() != 6*mesh.nodes.len() || row.iter().any(|v| !v.is_finite()))
        { return Err(bad("cavity trace requires every original physical shell mode")); }
        for selection in &surfaces.cavity {
            check_cancel(gate)?;
            let nodes = mesh.tris[selection.triangle];
            let (face, points, normal) = faces.remove(&key(nodes.map(|n| mesh.nodes[n])))
                .ok_or_else(|| bad("shell cavity facet is missing, duplicated or does not conform to the volume boundary"))?;
            let facet = mesh.facet(selection.triangle)
                .map_err(|_| bad("original shell cavity facet is invalid"))?;
            let sign = match selection.side { ShellFace::Positive => 1.0, ShellFace::Negative => -1.0 };
            let alignment: f64 = facet.frame[2].iter().zip(normal).map(|(a,b)| sign*a*b).sum();
            if alignment < 1.0 - 256.0*f64::EPSILON {
                return Err(bad("shell pressure ownership is not outward from the original cavity volume"));
            }
            let mut ordered = [0usize; 3];
            for (index, point) in points.iter().enumerate() {
                ordered[index] = *nodes.iter().find(|&&node| mesh.nodes[node] == *point)
                    .ok_or_else(|| bad("conforming cavity lost a shared physical vertex"))?;
            }
            for (trace, row) in traces.iter_mut().zip(rows) {
                let values = ordered.map(|node| (0..3).map(|axis| normal[axis]*row[6*node+axis]).sum());
                trace.push((face, values));
            }
        }
    }
    for opening in openings {
        for triangle in &opening.triangles {
            if faces.remove(&key(*triangle)).is_none() {
                return Err(bad("cavity mouth overlaps a solid pressure face or another aperture"));
            }
        }
    }
    if !faces.is_empty() || !area.is_finite() || area <= 0.0 {
        return Err(bad("every cavity boundary must have one original solid or explicit fluid-aperture owner"));
    }
    Ok((traces, area))
}

