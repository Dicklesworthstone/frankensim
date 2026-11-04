//! Source-authored longitudinal Maxwell memory in native frame AC and poles.
//!
//! The original member supplies E_infinity. Every positive layer arm contributes
//! its original delta_E and tau through the geometry-owned axial and two bending
//! strain traces. All native frame coordinates, memory states and optional
//! boundary states/constraints enter the SAME coupled pHS preparation.
//!
//! In the elastic-strain chart z=sqrt(k)*(extension-viscous_extension),
//! z_dot=-z/tau+sqrt(k)*extension_dot and opposing effort=sqrt(k)*z.
//! H=z*z/2 and heat=z*z/tau. Eliminating z only for interpretation gives
//! delta_K(s)=k*s/(s+1/tau); the implementation never eliminates these states.

mod memory;

use super::coupled::{
    FrameCoupledControls, FrameCoupledError, FrameDynamicPortLoad, FramePortDynamics,
    PreparedFrameResonance, limit,
};
use super::{Frame2d, Rayleigh};
use fs_exec::Cx;
use fs_material::section::relaxation::{ConcentricPronySection, LayerRelaxationProperties};
use fs_material::state_point::MaterialStatePointError;
use fs_math::c64::C64;
use fs_phs::harmonic::{HarmonicResponse, MAX_HARMONIC_DIMENSION};
use fs_phs::{DescriptorPortHamiltonian, QuadraticStorage};
use memory::{MemoryCoefficient, longitudinal_coefficients};
use std::collections::BTreeSet;

/// One original scalar longitudinal E(t) per bonded layer, explicitly governing
/// BOTH axial strain and bending. This is not a bending-only or torsional law.
/// Selectors resolve against the ORIGINAL material-bound frame member.
#[derive(Clone, Debug)]
pub struct FrameLongitudinalRelaxation {
    pub element: usize,
    pub layers: Vec<LayerRelaxationProperties>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameStrainChannel {
    Axial,
    SymmetricBending,
    AntisymmetricBending,
}
impl FrameStrainChannel {
    fn index(self) -> usize {
        match self {
            Self::Axial => 0,
            Self::SymmetricBending => 1,
            Self::AntisymmetricBending => 2,
        }
    }
}

/// Source receipts and all original arms, including exactly zero increments.
#[derive(Clone, Debug)]
pub struct FrameRelaxationMember {
    pub element: usize,
    pub source: ConcentricPronySection,
    dofs: [usize; 6],
    weights: [[f64; 6]; 3],
}

/// One positive source arm/channel. A zero projection onto free motion still
/// has its own energy and decay state; only an exactly zero increment is omitted.
#[derive(Clone, Copy, Debug)]
pub struct FrameRelaxationCoordinate {
    pub element: usize,
    pub layer: usize,
    pub source_arm: usize,
    pub channel: FrameStrainChannel,
    /// Index in the complete native `[q,p,material_memory,boundary]` state.
    pub native_state: usize,
    pub stiffness_n_m: f64,
    pub relaxation_time_s: f64,
    member: usize,
    coefficient: MemoryCoefficient,
}

#[derive(Debug)]
pub enum FrameRelaxationError {
    Source(MaterialStatePointError),
    Native(FrameCoupledError),
}
impl core::fmt::Display for FrameRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-FRAME-RELAXATION: {self:?}")
    }
}
impl std::error::Error for FrameRelaxationError {}
impl From<FrameCoupledError> for FrameRelaxationError {
    fn from(error: FrameCoupledError) -> Self {
        Self::Native(error)
    }
}
impl From<MaterialStatePointError> for FrameRelaxationError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}
fn invalid(message: &'static str) -> FrameRelaxationError {
    FrameCoupledError::Invalid(message).into()
}
fn poll(cx: &Cx<'_>) -> Result<(), FrameRelaxationError> {
    cx.checkpoint()
        .map_err(|_| FrameCoupledError::Cancelled.into())
}

/// Complete immutable native resonance with source identity and physical memory
/// observations. The same coupled owner performs factorization and pole search.
pub struct PreparedFrameRelaxation {
    coupled: PreparedFrameResonance,
    members: Vec<FrameRelaxationMember>,
    coordinates: Vec<FrameRelaxationCoordinate>,
}

#[derive(Clone, Debug)]
pub struct FrameRelaxationPhasor {
    pub coordinate: FrameRelaxationCoordinate,
    /// Original integrated strain coordinate [m], not point fibre strain.
    pub extension_m: C64,
    pub elastic_extension_m: C64,
    pub viscous_extension_m: C64,
    /// Opposing scalar force, mapped to nodal forces/moments by the SAME trace.
    pub opposing_force_n: C64,
    pub mean_stored_energy_j: f64,
    pub mean_dissipated_power_w: f64,
}

/// Only this constitutive mechanism's storage and heat. Native total response
/// also includes equilibrium, kinetic, Rayleigh/Kelvin and boundary accounts.
#[derive(Clone, Debug)]
pub struct FrameRelaxationResponse {
    pub coordinates: Vec<FrameRelaxationPhasor>,
    pub mean_stored_energy_j: f64,
    pub mean_dissipated_power_w: f64,
}

impl Frame2d {
    /// Prepare full native resonance with original layer longitudinal memory.
    ///
    /// The section must be the original material-bound, linear, reference-frame
    /// source. E_infinity is matched by the existing source resolver; it is not
    /// added a second time. Nonzero material Kelvin–Voigt loss on a selected
    /// member refuses to avoid counting an unspecified second calibration.
    /// Separately authored Rayleigh damping and unselected members remain active.
    ///
    /// Each positive layer arm supplies three physical coordinates. All original
    /// tau values and channel couplings survive, even if supports make a memory
    /// invisible to motion. Exact-zero arms retain source identity without an
    /// invented energy state. An optional boundary keeps all of its original
    /// states, mutual matrices, ports and trailing algebraic multipliers.
    ///
    /// State order is `[q,p,material_memory,original_boundary]`. q/p follow the
    /// native free-DOF order. Memory uses energy-normalized ELASTIC extension;
    /// it is not the viscous-extension chart used by transient relaxation APIs.
    ///
    /// Real-frequency AC is gated by the intersection of original source bands
    /// and native material-viscosity support. Pole queries describe the declared
    /// continuous Maxwell law's complex continuation; finite real-frequency
    /// source bands do not certify that continuation. The boundary source owns
    /// its additional applicability. Fixed temperature, reference density and
    /// geometry, small strain, perfect bonding and Euler–Bernoulli kinematics
    /// apply; no transient history, live thermal update or heat deposition occurs.
    ///
    /// Bounds include assembly of the constitutive/boundary block AND shared
    /// native preparation. Before resolving sources, three states per declared
    /// arm (including zero increments) are conservatively charged. Source
    /// metadata and caller-owned inputs retain the coupled controls' exclusion.
    pub fn prepare_relaxing_resonance(
        &self,
        rayleigh: Rayleigh,
        selections: &[FrameLongitudinalRelaxation],
        boundary: Option<FrameDynamicPortLoad<'_>>,
        controls: FrameCoupledControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedFrameRelaxation, FrameRelaxationError> {
        poll(cx)?;
        if selections.is_empty() || selections.len() > self.elements.len() {
            return Err(invalid(
                "nonempty unique original member selections required",
            ));
        }
        let n = self.ndof();
        let nf = self.fixed.iter().filter(|&&fixed| !fixed).count();
        limit(
            "full frame coordinates",
            n as u128,
            controls.maximum_coordinates as u128,
        )?;
        let external = boundary.as_ref().map(|load| load.dynamics.parts());
        let (nb, db, pb) = external.map_or((0, 0, 0), |(n, d, p, ..)| (n, d, p));
        let material_ports = selections
            .len()
            .checked_mul(3)
            .ok_or_else(|| invalid("material port count overflow"))?;
        let ports = material_ports
            .checked_add(pb)
            .ok_or_else(|| invalid("port count overflow"))?;
        limit(
            "material and boundary ports",
            ports as u128,
            controls.maximum_ports as u128,
        )?;
        let mut arms = 0usize;
        let mut unique = BTreeSet::new();
        for selection in selections {
            poll(cx)?;
            if !unique.insert(selection.element) {
                return Err(invalid("member relaxation was declared twice"));
            }
            for layer in &selection.layers {
                poll(cx)?;
                arms = arms
                    .checked_add(layer.branches.len())
                    .ok_or_else(|| invalid("source arm count overflow"))?;
            }
        }
        let memory_cap = 3_u128 * arms as u128;
        let load_cap = memory_cap + nb as u128;
        limit(
            "frame plus all declared memory and boundary states",
            2 * nf as u128 + load_cap,
            controls.maximum_states.min(MAX_HARMONIC_DIMENSION) as u128,
        )?;
        // Two PSD admissions plus block copies. Reserve this additional payload
        // and work before any source evaluation or dense block construction.
        let extra_bytes = 256 * load_cap.pow(2)
            + 64 * (load_cap + n as u128) * ports as u128
            + 1024 * (memory_cap + selections.len() as u128);
        let extra_work = 2048 * load_cap.pow(3)
            + 64 * load_cap.pow(2) * ports as u128
            + 4096 * (memory_cap + selections.len() as u128);
        limit(
            "constitutive composition bytes",
            extra_bytes,
            controls.maximum_payload_bytes as u128,
        )?;
        limit(
            "constitutive composition work",
            extra_work,
            controls.maximum_work as u128,
        )?;
        let native_controls = FrameCoupledControls {
            maximum_payload_bytes: controls.maximum_payload_bytes - extra_bytes as usize,
            maximum_work: controls.maximum_work - extra_work as u64,
            ..controls
        };
        // Validate a supplied boundary before indexing any original matrix.
        let external_q = if let Some((nb, db, pb, j, r, g, storage)) = external {
            if db > nb
                || j.len() != nb * nb
                || r.len() != nb * nb
                || g.len() != nb * pb
                || boundary
                    .is_none_or(|load| load.force_map.len() as u128 != n as u128 * pb as u128)
                || j.iter().chain(r).chain(g).any(|x| !x.is_finite())
                || boundary.is_some_and(|load| load.force_map.iter().any(|x| !x.is_finite()))
            {
                return Err(invalid(
                    "complete finite original boundary structure/map required",
                ));
            }
            let part = storage
                .quadratic_part()
                .ok_or_else(|| invalid("boundary must declare its complete quadratic storage"))?;
            if !part.exact
                || part.matrix.len() != nb * nb
                || part.matrix.iter().any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "complete finite quadratic boundary storage required",
                ));
            }
            for i in 0..nb {
                poll(cx)?;
                for j in 0..nb {
                    if (i >= db || j >= db) && part.matrix[i * nb + j] != 0.0 {
                        return Err(invalid(
                            "boundary algebraic multipliers cannot carry storage",
                        ));
                    }
                }
            }
            Some(part.matrix)
        } else {
            None
        };
        let mut members = Vec::with_capacity(selections.len());
        let mut coordinates = Vec::with_capacity(memory_cap as usize);
        let mut band = (0.0_f64, f64::MAX);
        for selection in selections {
            poll(cx)?;
            if self
                .material_viscous_section(selection.element)
                .is_some_and(|section| {
                    section.axial().viscous_rigidity_n_s() != 0.0
                        || section.bending().viscous_stiffness_n_m2_s() != 0.0
                })
            {
                return Err(invalid(
                    "Prony member already carries nonzero material Kelvin-Voigt loss",
                ));
            }
            let source = self
                .material_section(selection.element)
                .ok_or_else(|| {
                    invalid("relaxation requires the original material-bound frame section")
                })?
                .clone()
                .with_prony_bending(&selection.layers)?;
            let source_band = source.omega_band_rad_s();
            band = (band.0.max(source_band.0), band.1.min(source_band.1));
            if band.0 > band.1 {
                return Err(invalid("original source frequency supports do not overlap"));
            }
            let (dofs, weights, length) = self
                .reference_strain_traces(selection.element)
                .map_err(FrameCoupledError::Frame)?;
            for (layer, (geometry, law)) in source
                .elastic()
                .layer_integrals()
                .iter()
                .zip(source.spectra())
                .enumerate()
            {
                for (source_arm, &(delta, tau)) in law.terms.iter().enumerate() {
                    poll(cx)?;
                    if delta == 0.0 {
                        continue;
                    }
                    let coefficients = longitudinal_coefficients(
                        delta,
                        geometry.area_m2(),
                        geometry.second_moment_m4(),
                        length,
                        tau,
                    )
                    .ok_or_else(|| {
                        invalid("positive longitudinal memory coefficients are not representable")
                    })?;
                    for (channel, coefficient) in [
                        FrameStrainChannel::Axial,
                        FrameStrainChannel::SymmetricBending,
                        FrameStrainChannel::AntisymmetricBending,
                    ]
                    .into_iter()
                    .zip(coefficients)
                    {
                        coordinates.push(FrameRelaxationCoordinate {
                            element: selection.element,
                            layer,
                            source_arm,
                            channel,
                            native_state: 2 * nf + coordinates.len(),
                            stiffness_n_m: coefficient.stiffness,
                            relaxation_time_s: tau,
                            member: members.len(),
                            coefficient,
                        });
                    }
                }
            }
            members.push(FrameRelaxationMember {
                element: selection.element,
                source,
                dofs,
                weights,
            });
        }
        let nm = coordinates.len();
        let nl = nm + nb;
        let dl = nm + db;
        let mut q = vec![0.0; dl * dl];
        let mut j = vec![0.0; nl * nl];
        let mut r = vec![0.0; nl * nl];
        let mut g = vec![0.0; nl * ports];
        let mut force_map = vec![0.0; n * ports];
        for (member, source) in members.iter().enumerate() {
            poll(cx)?;
            for (channel, weights) in source.weights.iter().enumerate() {
                for (&dof, &weight) in source.dofs.iter().zip(weights) {
                    force_map[dof * ports + 3 * member + channel] = weight;
                }
            }
        }
        for (i, coordinate) in coordinates.iter().enumerate() {
            poll(cx)?;
            q[i * dl + i] = 1.0;
            r[i * nl + i] = coordinate.coefficient.rate;
            g[i * ports + 3 * coordinate.member + coordinate.channel.index()] =
                coordinate.coefficient.root;
        }
        if let Some((_, _, _, jb, rb, gb, _)) = external {
            let qb = external_q
                .as_ref()
                .expect("validated original boundary storage");
            for i in 0..nb {
                poll(cx)?;
                for k in 0..nb {
                    j[(nm + i) * nl + nm + k] = jb[i * nb + k];
                    r[(nm + i) * nl + nm + k] = rb[i * nb + k];
                    if i < db && k < db {
                        q[(nm + i) * dl + nm + k] = qb[i * nb + k];
                    }
                }
                g[(nm + i) * ports + material_ports..(nm + i + 1) * ports]
                    .copy_from_slice(&gb[i * pb..(i + 1) * pb]);
            }
            let map = boundary.expect("validated original boundary").force_map;
            for i in 0..n {
                poll(cx)?;
                force_map[i * ports + material_ports..(i + 1) * ports]
                    .copy_from_slice(&map[i * pb..(i + 1) * pb]);
            }
        }
        poll(cx)?;
        let storage = QuadraticStorage::new(q, dl).map_err(FrameCoupledError::Physical)?;
        poll(cx)?;
        let load = DescriptorPortHamiltonian::new(nl, dl, ports, j, r, g, Box::new(storage))
            .map_err(FrameCoupledError::Physical)?;
        poll(cx)?;
        let mut coupled = self.prepare_coupled_resonance(
            rayleigh,
            FrameDynamicPortLoad {
                force_map: &force_map,
                dynamics: FramePortDynamics::Descriptor(&load),
            },
            native_controls,
            cx,
        )?;
        coupled.intersect_material_band(band)?;
        poll(cx)?;
        Ok(PreparedFrameRelaxation {
            coupled,
            members,
            coordinates,
        })
    }
}

impl PreparedFrameRelaxation {
    /// Same complete coupled AC/local-pole/contour owner; the original material
    /// frequency gate is already installed and cannot be bypassed through it.
    pub fn coupled(&self) -> &PreparedFrameResonance {
        &self.coupled
    }
    pub fn members(&self) -> &[FrameRelaxationMember] {
        &self.members
    }
    pub fn coordinates(&self) -> &[FrameRelaxationCoordinate] {
        &self.coordinates
    }
    pub fn memory_state_range(&self) -> core::ops::Range<usize> {
        let start = 2 * self.coupled.free_dofs().len();
        start..start + self.coordinates.len()
    }
    /// Every original boundary state in its unchanged source order; empty when
    /// no boundary is supplied. Its algebraic range is the coupled owner's range.
    pub fn boundary_state_range(&self) -> core::ops::Range<usize> {
        self.memory_state_range().end..self.coupled.dimension()
    }

    /// Observe peak phasors from THIS model's native forced response. This
    /// reconstructs original extension, elastic/viscous partition, storage and
    /// heat without solving a second constitutive response. Returned component
    /// heat is already included in native total dissipation; do not add it twice.
    /// Dimensional/source-band checks do not authenticate a foreign response.
    pub fn harmonic_memory(
        &self,
        response: &HarmonicResponse,
        cx: &Cx<'_>,
    ) -> Result<FrameRelaxationResponse, FrameRelaxationError> {
        poll(cx)?;
        let band = self
            .coupled
            .material_omega_band_rad_s()
            .expect("source material band");
        if response.state.len() != self.coupled.dimension()
            || response
                .state
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
            || !response.omega_rad_s.is_finite()
            || response.omega_rad_s < band.0
            || response.omega_rad_s > band.1
        {
            return Err(invalid(
                "finite original state and source-supported forcing frequency required",
            ));
        }
        let mut extensions = Vec::with_capacity(self.members.len());
        for member in &self.members {
            poll(cx)?;
            let mut extension = [C64::ZERO; 3];
            for (channel, weights) in member.weights.iter().enumerate() {
                for (&dof, &weight) in member.dofs.iter().zip(weights) {
                    if let Ok(index) = self.coupled.free_dofs().binary_search(&dof) {
                        extension[channel] =
                            extension[channel] + response.state[index].scale(weight);
                    }
                }
            }
            extensions.push(extension);
        }
        let mut result = FrameRelaxationResponse {
            coordinates: Vec::with_capacity(self.coordinates.len()),
            mean_stored_energy_j: 0.0,
            mean_dissipated_power_w: 0.0,
        };
        for &coordinate in &self.coordinates {
            poll(cx)?;
            let extension = extensions[coordinate.member][coordinate.channel.index()];
            let seen = coordinate
                .coefficient
                .observe(extension, response.state[coordinate.native_state])
                .ok_or_else(|| invalid("original memory observation is not representable"))?;
            result.mean_stored_energy_j += seen.energy;
            result.mean_dissipated_power_w += seen.heat;
            result.coordinates.push(FrameRelaxationPhasor {
                coordinate,
                extension_m: extension,
                elastic_extension_m: seen.elastic,
                viscous_extension_m: seen.viscous,
                opposing_force_n: seen.force,
                mean_stored_energy_j: seen.energy,
                mean_dissipated_power_w: seen.heat,
            });
        }
        if !result.mean_stored_energy_j.is_finite() || !result.mean_dissipated_power_w.is_finite() {
            return Err(invalid(
                "total material memory observation is not representable",
            ));
        }
        poll(cx)?;
        Ok(result)
    }
}
