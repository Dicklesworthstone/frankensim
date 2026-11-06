//! Complete structural displacement and heterogeneous pressure fields in one solve.
//!
//! The pressure potential satisfies `p = psi_dot`. With `B[a,s] =
//! integral(phi_a * normal_displacement_s dA)`, the original equations are
//!
//! `Ms u_ddot + Cs u_dot + Ks u - B^T p = F`,
//! `Ma psi_ddot + Cwall psi_dot + Ka psi + B u_dot = s`.
//!
//! `Ma = integral(phi_i phi_j / (rho c^2))` and
//! `Ka = integral(grad(phi_i).grad(phi_j) / rho)` are copied from the original
//! heterogeneous owner. The two interface blocks are exactly opposite
//! transposes, so their instantaneous work cancels. The canonical
//! `fs-modal` harmonic owner factors the whole pencil, including its skew
//! velocity coupling; no bare-fluid inverse or pressure-mode truncation is used.
//!
//! Every phasor here uses `exp(-i omega t)`, including prescribed flow. This
//! differs from the standalone heterogeneous pressure harmonic API. Geometry,
//! coefficients and real trace quadrature carry no phasor convention to convert.
//! All original pressure coordinates, including constant compression, remain.
//!
//! Scope: frozen linear fluids, small structural motion, pressure/normal-velocity
//! coupling and optional real nonnegative local admittance relative to the moving
//! wall. The solid adapter is conforming P1 on planar faces. This does not supply
//! tangential viscous traction, moving geometry, exterior radiation or a causal
//! realization of frequency-dependent wall data. Full spatial coordinates still
//! require mesh/order convergence; dense harmonic preparation is explicitly bounded.

/// Full original structural and pressure fields on one exact held-input timeline.
pub mod time;

use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::AcousticFeError;
use fs_feec::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_feec::acoustics::helmholtz::heterogeneous::{HeterogeneousHelmholtzFe, VelocityTrace};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::Csr;
use std::collections::{BTreeMap, BTreeSet};

/// A P1 outward-from-fluid normal displacement shape at the three face vertices.
/// Values follow the supplied vertex order, per unit original structural coordinate.
pub type NormalDisplacementTrace = ([u32; 3], [f64; 3]);

/// The original complete free-coordinate structural pencil, with no normalization.
/// A missing velocity operator means exactly zero damping/gyroscopic coupling.
pub struct StructuralMatrices<'s> {
    pub mass: &'s Csr,
    pub stiffness: &'s Csr,
    pub damping: Option<&'s Csr>,
}

/// One conforming fluid/solid face. Corresponding vertices must coincide exactly;
/// winding does not set the normal, which comes from the owning fluid volume.
#[derive(Clone, Copy, Debug)]
pub struct ConformingTetFace {
    pub fluid: [u32; 3],
    pub solid: [usize; 3],
}

/// Original solid assembly and the geometry from which it was assembled.
/// `free_dofs` supplies the actual constrained-to-free map; fixed coordinates
/// are never recreated as independent moving-wall coordinates.
pub struct ConformingTetStructure<'s> {
    pub assembly: &'s TetElasticAssembly,
    pub nodes_m: &'s [[f64; 3]],
    pub tetrahedra: &'s [[usize; 4]],
    pub interface: &'s [ConformingTetFace],
}

/// Original geometry and harmonic-owner refusals remain visible.
#[derive(Debug)]
pub enum StructuralInteriorError {
    Invalid(&'static str),
    Interior(AcousticFeError),
    Harmonic(HarmonicResponseError),
}

impl core::fmt::Display for StructuralInteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "structural interior: {what}"),
            Self::Interior(error) => error.fmt(f),
            Self::Harmonic(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for StructuralInteriorError {}

impl From<AcousticFeError> for StructuralInteriorError {
    fn from(error: AcousticFeError) -> Self {
        Self::Interior(error)
    }
}
impl From<HarmonicResponseError> for StructuralInteriorError {
    fn from(error: HarmonicResponseError) -> Self {
        Self::Harmonic(error)
    }
}

/// Complete original spatial pencil. Preparation reuses it across frequencies;
/// each prepared frequency reuses the canonical factor across independent drives.
pub struct StructuralInterior<'a, 'f, 'c> {
    interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
    structural_dofs: usize,
    mass: Vec<f64>,
    damping: Vec<f64>,
    stiffness: Vec<f64>,
    coupling: Vec<f64>,
    controls: HarmonicResponseControls,
}

impl<'a, 'f, 'c> StructuralInterior<'a, 'f, 'c> {
    /// Integrate one actual normal trace per original structural coordinate.
    ///
    /// All source matrix entries survive, including off-diagonal inertia and
    /// structural gyroscopic coupling. The harmonic owner admits SPD mass,
    /// PSD stiffness and the symmetric dissipative part at preparation.
    /// Unmoving coordinates use an empty trace. Resistive wall groups must have
    /// distinct faces; their admittance is relative to the structural motion.
    ///
    /// The coordinate cap covers structure plus pressure. The payload cap admits
    /// `176*n*n + 256*n` bytes for this dense pencil, its canonical preparation
    /// and one response, excluding caller sources, geometry and allocator overhead.
    pub fn new(
        interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
        structure: StructuralMatrices<'_>,
        normal_traces: &[Vec<NormalDisplacementTrace>],
        walls: &[ResistiveBoundary],
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<Self, StructuralInteriorError> {
        poll(cx)?;
        let ns = structure.mass.nrows();
        let na = interior.ndof();
        let n = admit(ns, na, controls)?;
        if normal_traces.len() != ns {
            return Err(invalid(
                "one normal trace collection per free structural coordinate required",
            ));
        }
        let mut mass = vec![0.0; n * n];
        let mut stiffness = vec![0.0; n * n];
        let mut damping = vec![0.0; n * n];
        copy_block(&mut mass, n, 0, structure.mass, ns, cx)?;
        copy_block(&mut stiffness, n, 0, structure.stiffness, ns, cx)?;
        if let Some(source) = structure.damping {
            copy_block(&mut damping, n, 0, source, ns, cx)?;
        }
        copy_block(&mut mass, n, ns, interior.mass(), na, cx)?;
        copy_block(&mut stiffness, n, ns, interior.stiffness(), na, cx)?;
        let mut coupling = vec![0.0; na * ns];
        for (s, traces) in normal_traces.iter().enumerate() {
            poll(cx)?;
            let load = interior.geometry().boundary_linear_load(traces)?;
            for (a, &value) in load.iter().enumerate() {
                coupling[a * ns + s] = value;
                damping[(ns + a) * n + s] = value;
                damping[s * n + ns + a] = -value;
            }
        }
        let mut seen = BTreeSet::new();
        for wall in walls {
            poll(cx)?;
            if !wall.admittance_m_pa_s.is_finite() || wall.admittance_m_pa_s < 0.0 {
                return Err(invalid("wall admittance must be finite and nonnegative"));
            }
            for &face in &wall.faces {
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(invalid(
                        "resistive wall faces must be unique across all groups",
                    ));
                }
            }
            let surface = interior.geometry().boundary_mass(&wall.faces)?;
            for a in 0..na {
                let (cols, values) = surface.row(a);
                for (&b, &value) in cols.iter().zip(values) {
                    let target = &mut damping[(ns + a) * n + ns + b];
                    *target += wall.admittance_m_pa_s * value;
                    if !target.is_finite() {
                        return Err(invalid("wall admittance assembly overflow"));
                    }
                }
            }
        }
        poll(cx)?;
        Ok(Self {
            interior,
            structural_dofs: ns,
            mass,
            damping,
            stiffness,
            coupling,
            controls,
        })
    }

    /// Bind the actual P1 solid shape functions to the original pressure trace.
    ///
    /// Interface faces must belong to exactly one supplied solid tetrahedron,
    /// on the opposite side of the fluid. The original consistent mass,
    /// anisotropic stiffness and free-coordinate order come from `assembly`.
    /// The supplied geometry must be the source of that assembly.
    pub fn from_tet_elastic(
        interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
        solid: ConformingTetStructure<'_>,
        damping: Option<&Csr>,
        walls: &[ResistiveBoundary],
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<Self, StructuralInteriorError> {
        poll(cx)?;
        let ns = solid.assembly.free_dofs.len();
        admit(ns, interior.ndof(), controls)?;
        let full = solid
            .nodes_m
            .len()
            .checked_mul(3)
            .ok_or_else(|| invalid("solid coordinate extent overflow"))?;
        let mut free = BTreeMap::new();
        for (reduced, &original) in solid.assembly.free_dofs.iter().enumerate() {
            if original >= full || free.insert(original, reduced).is_some() {
                return Err(invalid(
                    "solid free-coordinate map must be unique and in range",
                ));
            }
        }
        let mut owners = BTreeMap::new();
        let mut fluid_faces = BTreeSet::new();
        for (index, face) in solid.interface.iter().enumerate() {
            let mut key = face.solid;
            key.sort_unstable();
            let mut fluid = face.fluid;
            fluid.sort_unstable();
            if key[0] == key[1]
                || key[1] == key[2]
                || key[2] >= solid.nodes_m.len()
                || owners.insert(key, (index, 0usize, 0usize)).is_some()
                || !fluid_faces.insert(fluid)
            {
                return Err(invalid(
                    "conforming interface requires distinct valid fluid and solid faces",
                ));
            }
        }
        for tet in solid.tetrahedra {
            poll(cx)?;
            let mut key = *tet;
            key.sort_unstable();
            if key[3] >= solid.nodes_m.len() || key.windows(2).any(|pair| pair[0] == pair[1]) {
                return Err(invalid("solid tetrahedron connectivity is invalid"));
            }
            for opposite in 0..4 {
                let mut face = [0; 3];
                let mut next = 0;
                for (local, &vertex) in tet.iter().enumerate() {
                    if local != opposite {
                        face[next] = vertex;
                        next += 1;
                    }
                }
                face.sort_unstable();
                if let Some((_, count, vertex)) = owners.get_mut(&face) {
                    *count += 1;
                    *vertex = tet[opposite];
                }
            }
        }
        let mut traces = vec![Vec::new(); ns];
        for (_, (index, count, opposite)) in owners {
            poll(cx)?;
            if count != 1 {
                return Err(invalid(
                    "each interface face must have exactly one solid volume owner",
                ));
            }
            let face = solid.interface[index];
            let (points, normal, _) = interior.geometry().boundary_geometry(face.fluid)?;
            for local in 0..3 {
                if points[local] != solid.nodes_m[face.solid[local]] {
                    return Err(invalid(
                        "fluid and solid interface vertices must coincide exactly",
                    ));
                }
            }
            let outside: f64 = (0..3)
                .map(|axis| normal[axis] * (solid.nodes_m[opposite][axis] - points[0][axis]))
                .sum();
            if !outside.is_finite() || outside <= 0.0 {
                return Err(invalid(
                    "solid interface owner must lie outside the fluid volume",
                ));
            }
            for axis in 0..3 {
                for local in 0..3 {
                    if let Some(&coordinate) = free.get(&(3 * face.solid[local] + axis)) {
                        let mut values = [0.0; 3];
                        values[local] = normal[axis];
                        traces[coordinate].push((face.fluid, values));
                    }
                }
            }
        }
        Self::new(
            interior,
            StructuralMatrices {
                mass: &solid.assembly.mass,
                stiffness: &solid.assembly.stiffness,
                damping,
            },
            &traces,
            walls,
            controls,
            cx,
        )
    }

    #[must_use]
    pub const fn structural_dofs(&self) -> usize {
        self.structural_dofs
    }

    #[must_use]
    pub fn pressure_dofs(&self) -> usize {
        self.interior.ndof()
    }

    /// Full row-major `pressure_dofs x structural_dofs` reciprocal work map.
    #[must_use]
    pub fn coupling(&self) -> &[f64] {
        &self.coupling
    }

    /// Project prescribed INWARD velocity, using the original full pressure
    /// trace quadrature. The result is the generalized volume-flow load [m^3/s]
    /// for `solve`, without area normalization or an extra density/frequency.
    pub fn project_inward_velocity(
        &self,
        traces: &[VelocityTrace],
    ) -> Result<Vec<C64>, StructuralInteriorError> {
        let real: Vec<_> = traces
            .iter()
            .map(|&(face, v)| (face, v.map(|z| z.re)))
            .collect();
        let imag: Vec<_> = traces
            .iter()
            .map(|&(face, v)| (face, v.map(|z| z.im)))
            .collect();
        let real = self.interior.geometry().boundary_linear_load(&real)?;
        let imag = self.interior.geometry().boundary_linear_load(&imag)?;
        Ok(real
            .into_iter()
            .zip(imag)
            .map(|(re, im)| C64::new(re, im))
            .collect())
    }

    /// Prepare the complete coupled system at one strictly positive frequency.
    /// A bare acoustic pole need not be invertible: the coupled system is factored
    /// directly. True lossless coupled poles retain the canonical refusal.
    pub fn prepare(
        &self,
        frequency_hz: f64,
        cx: &Cx<'_>,
    ) -> Result<PreparedStructuralInterior<'_, 'a, 'f, 'c>, StructuralInteriorError> {
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
            return Err(invalid(
                "strictly positive finite harmonic frequency required",
            ));
        }
        let harmonic = PreparedHarmonicResponse::try_new(
            &self.mass,
            &self.damping,
            &self.stiffness,
            self.structural_dofs + self.pressure_dofs(),
            core::f64::consts::TAU * frequency_hz,
            self.controls,
            cx,
        )?;
        Ok(PreparedStructuralInterior {
            model: self,
            harmonic,
        })
    }
}

/// The original whole-system factor at one frequency.
pub struct PreparedStructuralInterior<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    harmonic: PreparedHarmonicResponse,
}

/// All fields remain in their original free/pressure coordinate order.
#[derive(Clone, Debug)]
pub struct StructuralInteriorResponse {
    pub structural_displacement: Vec<C64>,
    pub structural_velocity: Vec<C64>,
    pub pressure_potential_pa_s: Vec<C64>,
    pub pressure_pa: Vec<C64>,
    /// Generalized force exerted by the fluid on the structure, `B^T p`.
    pub fluid_force: Vec<C64>,
    /// Original pressure-basis integrals of structural outward normal velocity.
    pub displaced_flow_m3_s: Vec<C64>,
    pub average_structural_input_power_w: f64,
    pub average_acoustic_input_power_w: f64,
    /// Signed interface work, positive from fluid into structure.
    pub average_fluid_to_structure_power_w: f64,
    /// Complete symmetric structural velocity-operator work.
    pub average_structural_dissipation_w: f64,
    /// Independently integrated relative-wall admittance work.
    pub average_wall_dissipation_w: f64,
    /// Whole-system input minus dissipation, from the canonical harmonic owner.
    pub power_balance_residual_w: f64,
    pub relative_residual: f64,
}

impl PreparedStructuralInterior<'_, '_, '_, '_> {
    /// Apply peak original structural forces and projected inward volume flow.
    /// The flow load is `s_a = integral(phi_a v_in dA)` [m^3/s], or an original
    /// projected volume-injection source in the same basis. Its conjugate
    /// pressure is the returned `pressure_pa`, never the potential.
    pub fn solve(
        &self,
        structural_force: &[C64],
        inward_flow_load: &[C64],
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorResponse, StructuralInteriorError> {
        poll(cx)?;
        let ns = self.model.structural_dofs;
        let na = self.model.pressure_dofs();
        if structural_force.len() != ns || inward_flow_load.len() != na {
            return Err(invalid(
                "force and flow loads must cover every original coordinate",
            ));
        }
        let mut force = Vec::with_capacity(ns + na);
        force.extend_from_slice(structural_force);
        force.extend_from_slice(inward_flow_load);
        let mut response = self.harmonic.solve(&force, cx)?;
        let potential = response.displacement.split_off(ns);
        let pressure = response.velocity.split_off(ns);
        let mut fluid_force = vec![C64::ZERO; ns];
        let mut flow = vec![C64::ZERO; na];
        for a in 0..na {
            poll(cx)?;
            for s in 0..ns {
                let b = self.model.coupling[a * ns + s];
                fluid_force[s] = fluid_force[s] + pressure[a].scale(b);
                flow[a] = flow[a] + response.velocity[s].scale(b);
            }
        }
        let structural_power = work(structural_force, &response.velocity);
        let acoustic_power = work(inward_flow_load, &pressure);
        let exchange = work(&fluid_force, &response.velocity);
        let structural_loss = dissipation(&self.model.damping, ns + na, 0, &response.velocity);
        let wall_loss = dissipation(&self.model.damping, ns + na, ns, &pressure);
        if fluid_force
            .iter()
            .chain(&flow)
            .any(|z| !z.re.is_finite() || !z.im.is_finite())
            || [
                structural_power,
                acoustic_power,
                exchange,
                structural_loss,
                wall_loss,
            ]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(invalid("interface force, flow or work overflow"));
        }
        Ok(StructuralInteriorResponse {
            structural_displacement: response.displacement,
            structural_velocity: response.velocity,
            pressure_potential_pa_s: potential,
            pressure_pa: pressure,
            fluid_force,
            displaced_flow_m3_s: flow,
            average_structural_input_power_w: structural_power,
            average_acoustic_input_power_w: acoustic_power,
            average_fluid_to_structure_power_w: exchange,
            average_structural_dissipation_w: structural_loss,
            average_wall_dissipation_w: wall_loss,
            power_balance_residual_w: response.power_balance_residual_w,
            relative_residual: response.relative_residual,
        })
    }
}

fn invalid(what: &'static str) -> StructuralInteriorError {
    StructuralInteriorError::Invalid(what)
}

fn poll(cx: &Cx<'_>) -> Result<(), StructuralInteriorError> {
    cx.checkpoint()
        .map_err(|_| HarmonicResponseError::Cancelled.into())
}

fn admit(
    ns: usize,
    na: usize,
    controls: HarmonicResponseControls,
) -> Result<usize, StructuralInteriorError> {
    let n = ns
        .checked_add(na)
        .ok_or_else(|| invalid("coupled coordinate extent overflow"))?;
    if ns == 0 || na == 0 {
        return Err(invalid(
            "nonempty original structural and pressure spaces required",
        ));
    }
    if n > controls.maximum_coordinates {
        return Err(HarmonicResponseError::Budget {
            quantity: "structural-interior coordinates",
            requested: n,
            cap: controls.maximum_coordinates,
        }
        .into());
    }
    let bytes = n
        .checked_mul(n)
        .and_then(|v| v.checked_mul(176))
        .and_then(|v| n.checked_mul(256).and_then(|w| v.checked_add(w)))
        .ok_or_else(|| invalid("coupled payload extent overflow"))?;
    if bytes > controls.maximum_payload_bytes {
        return Err(HarmonicResponseError::Budget {
            quantity: "structural-interior numerical payload bytes",
            requested: bytes,
            cap: controls.maximum_payload_bytes,
        }
        .into());
    }
    Ok(n)
}

fn copy_block(
    target: &mut [f64],
    n: usize,
    offset: usize,
    source: &Csr,
    size: usize,
    cx: &Cx<'_>,
) -> Result<(), StructuralInteriorError> {
    if source.nrows() != size || source.ncols() != size {
        return Err(invalid(
            "source matrices must use one complete original coordinate space",
        ));
    }
    for row in 0..size {
        poll(cx)?;
        let (cols, values) = source.row(row);
        for (&col, &value) in cols.iter().zip(values) {
            if !value.is_finite() {
                return Err(invalid("source matrix coefficients must be finite"));
            }
            target[(offset + row) * n + offset + col] = value;
        }
    }
    Ok(())
}

fn work(effort: &[C64], flow: &[C64]) -> f64 {
    0.5 * effort
        .iter()
        .zip(flow)
        .map(|(e, f)| e.re * f.re + e.im * f.im)
        .sum::<f64>()
}

fn dissipation(damping: &[f64], n: usize, offset: usize, velocity: &[C64]) -> f64 {
    let mut loss = 0.0;
    for i in 0..velocity.len() {
        for j in 0..velocity.len() {
            let a = damping[(offset + i) * n + offset + j];
            let b = damping[(offset + j) * n + offset + i];
            loss += (0.5 * a + 0.5 * b)
                * (velocity[i].re * velocity[j].re + velocity[i].im * velocity[j].im);
        }
    }
    0.5 * loss
}

#[cfg(test)]
mod tests;
