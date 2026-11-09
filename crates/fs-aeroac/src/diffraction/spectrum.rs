//! World-space, coherent illumination of the original exact straight edge.
//!
//! The Helmholtz equation and either ideal boundary condition are linear, so
//! a supplied propagating angular spectrum is scattered by summing the exact
//! Sommerfeld fields, including their velocity, before forming any intensity.
//! The caller owns spectral quadrature and its convergence, not this operator.
//! For the incident angular-spectrum convention see §2.1 of
//! <https://pmc.ncbi.nlm.nih.gov/articles/PMC7567611/>. No paraxial expansion,
//! angular interpolation, ray visibility test or order truncation is inserted.

use super::{HalfPlaneDiffraction, cis, finite, invalid, poll};
use crate::AeroacError;
use crate::curle3d::{CoherentCompactField, IdealPlaneBoundary};
use crate::observation::HarmonicFieldLine;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

/// Hard storage/work cap, in addition to the caller's component budget.
pub const MAX_HALF_PLANE_COMPONENTS: usize = 4096;

const SCOPE: &str = "Prescribed coherent propagating angular spectrum diffracted by one \
placed infinite straight-edged, zero-thickness rigid or pressure-release half-plane in a \
stationary homogeneous lossless fluid. Each supplied component uses exact Sommerfeld \
pressure and Euler velocity; the caller owns angular quadrature convergence. No evanescent \
incidence, finite-source reconstruction, finite barriers, wall impedance, general wedges, \
flow, absorption, source feedback or measured-SPL qualification.";

const FLOW_SCOPE: &str = "Prescribed coherent propagating angular spectrum diffracted by \
one placed ideal half-plane with uniform subsonic flow parallel to its edge. Each \
direction retains its own intrinsic frequency and wave number on one laboratory clock; \
coherent pressure and perturbation velocity are summed before pressure work or advected \
energy. No cross-edge flow, shear, finite barrier, wall impedance, absorption, source \
feedback or measured-SPL qualification.";

/// A placed half-plane: origin on its edge, z along the edge and y its normal.
///
/// The sheet extends in the positive x direction, `x = y cross z`. Both fluid
/// faces are retained. Inputs need finite nonzero directions orthogonal within
/// 64 epsilons after normalization; larger nonorthogonality refuses. The stored
/// right-handed basis removes only admitted normalization/roundoff error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HalfPlaneFrame {
    origin_m: [f64; 3],
    axes: [[f64; 3]; 3],
}

impl HalfPlaneFrame {
    /// Bind actual world geometry without inferring a boundary or material.
    pub fn new(
        origin_m: [f64; 3],
        edge_direction: [f64; 3],
        positive_sheet_normal: [f64; 3],
    ) -> Result<Self, AeroacError> {
        if origin_m.iter().any(|v| !v.is_finite()) {
            return Err(invalid("half-plane frame needs a finite origin"));
        }
        let z = normalized(edge_direction)?;
        let normal = normalized(positive_sheet_normal)?;
        if dot(z, normal).to_f64().abs() > 64.0 * f64::EPSILON {
            return Err(invalid(
                "half-plane edge and sheet normal must be orthogonal",
            ));
        }
        let x = normalized(cross(normal, z))?;
        let y = normalized(cross(z, x))?;
        Ok(Self {
            origin_m,
            axes: [x, y, z],
        })
    }

    pub const fn origin_m(&self) -> [f64; 3] {
        self.origin_m
    }

    /// Unit world vectors of the local x, y and z axes, in that order.
    pub const fn axes(&self) -> [[f64; 3]; 3] {
        self.axes
    }

    fn local_direction(&self, world: [f64; 3]) -> [f64; 3] {
        self.axes.map(|axis| dot(axis, world).to_f64())
    }

    pub(super) fn local_receiver(&self, receiver: HalfPlaneReceiver) -> Result<[f64; 3], AeroacError> {
        match receiver {
            HalfPlaneReceiver::WorldPosition(point) => {
                if point.iter().any(|v| !v.is_finite()) {
                    return Err(invalid(
                        "half-plane receiver needs finite world coordinates",
                    ));
                }
                let local = self
                    .axes
                    .map(|axis| translated_dot(axis, point, self.origin_m).to_f64());
                if local[0] >= 0.0 && local[1] == 0.0 {
                    return Err(invalid("on-sheet receiver must explicitly select its face"));
                }
                Ok(local)
            }
            HalfPlaneReceiver::SheetFace {
                distance_from_edge_m,
                along_edge_m,
                face,
            } => {
                if !distance_from_edge_m.is_finite()
                    || distance_from_edge_m <= 0.0
                    || !along_edge_m.is_finite()
                {
                    return Err(invalid("half-plane face receiver must lie beyond its edge"));
                }
                let y = match face {
                    HalfPlaneFace::PositiveNormal => 0.0,
                    HalfPlaneFace::NegativeNormal => -0.0,
                };
                Ok([distance_from_edge_m, y, along_edge_m])
            }
        }
    }
}

/// Select one of the two physically distinct limiting sheet faces.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HalfPlaneFace {
    PositiveNormal,
    NegativeNormal,
}

/// Receiver geometry, with an explicit representation for exact surface traces.
///
/// World-coordinate rounding cannot reliably encode the sign of zero after a
/// rotation. Surface traces therefore use physical sheet coordinates and a
/// named face. World points that resolve exactly onto the sheet refuse; no
/// tolerance snaps a nearby fluid point onto it.
#[derive(Clone, Copy, Debug)]
pub enum HalfPlaneReceiver {
    WorldPosition([f64; 3]),
    SheetFace {
        /// Positive distance along the sheet away from the edge [m].
        distance_from_edge_m: f64,
        /// Signed displacement along the edge from the frame origin [m].
        along_edge_m: f64,
        face: HalfPlaneFace,
    },
}

/// One coherent propagating component at the common world reference point.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneIncidentWave {
    /// Finite unit direction of travel in world coordinates. Its resolved
    /// component along the sheet normal must exceed 64 epsilons in magnitude.
    pub propagation_direction: [f64; 3],
    /// Complex peak pressure [Pa], exp(-i omega t), at the declared reference.
    /// Include any angular quadrature weight here; no normalization is added.
    pub pressure_pa: C64,
}

/// Explicit bound on preparation storage and work per receiver.
#[derive(Clone, Copy, Debug)]
pub struct HalfPlaneSpectrumBudget {
    /// Maximum supplied components, also limited by MAX_HALF_PLANE_COMPONENTS.
    pub max_components: usize,
}

/// Reusable exact-edge operator for one material, frequency and coherent spectrum.
///
/// Preparation stores O(N) original single-wave kernels and phases amplitudes
/// from the explicit world reference to the edge origin exactly once. Receiver
/// evaluation allocates nothing, costs O(N), and retains every supplied mode.
/// Cancellation is checked per component and within each Fresnel evaluation.
/// This is exact for the supplied finite propagating spectrum; quadrature of a
/// continuous spectrum requires a caller-measured angular refinement ladder.
#[derive(Clone, Debug)]
pub struct PreparedHalfPlaneSpectrum {
    frame: HalfPlaneFrame,
    reference_point_m: [f64; 3],
    components: Vec<(HalfPlaneDiffraction, C64)>,
}

impl PreparedHalfPlaneSpectrum {
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        frame: HalfPlaneFrame,
        boundary: IdealPlaneBoundary,
        reference_point_m: [f64; 3],
        incident: &[HalfPlaneIncidentWave],
        budget: HalfPlaneSpectrumBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::with_edge_parallel_flow(
            acoustic,
            frequency_hz,
            frame,
            boundary,
            0.0,
            reference_point_m,
            incident,
            budget,
            cancelled,
        )
    }

    /// Prepare one laboratory-frequency spectrum in a uniform flow along the
    /// placed edge. Speed is signed along `frame.axes()[2]`; the incident
    /// directions are wave normals, not laboratory group-velocity directions.
    /// Each component has its own `k = omega/(c + U*n_edge)`, including in
    /// the reference-point phase. A common scalar wave number would erase
    /// physical interference between upstream and downstream components.
    pub fn with_edge_parallel_flow(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        frame: HalfPlaneFrame,
        boundary: IdealPlaneBoundary,
        edge_flow_velocity_m_s: f64,
        reference_point_m: [f64; 3],
        incident: &[HalfPlaneIncidentWave],
        budget: HalfPlaneSpectrumBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if incident.is_empty()
            || incident.len() > budget.max_components
            || incident.len() > MAX_HALF_PLANE_COMPONENTS
        {
            return Err(invalid(
                "half-plane angular spectrum is empty or exceeds its component budget",
            ));
        }
        if reference_point_m.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "half-plane incident spectrum needs a finite phase reference",
            ));
        }
        let mut components = Vec::new();
        components
            .try_reserve_exact(incident.len())
            .map_err(|_| invalid("half-plane angular spectrum storage could not be reserved"))?;
        for component in incident {
            poll(&mut cancelled)?;
            if !finite(component.pressure_pa) {
                return Err(invalid("half-plane incident pressure must be finite"));
            }
            let direction = frame.local_direction(component.propagation_direction);
            if direction[1].abs() <= 64.0 * f64::EPSILON {
                return Err(invalid(
                    "half-plane incidence is grazing or unresolved in the world frame",
                ));
            }
            let kernel = HalfPlaneDiffraction::with_edge_parallel_flow(
                acoustic,
                frequency_hz,
                direction,
                boundary,
                edge_flow_velocity_m_s,
            )?;
            let phase = Dd::from_f64(kernel.k)
                * translated_dot(
                    component.propagation_direction,
                    frame.origin_m,
                    reference_point_m,
                );
            if !phase.to_f64().is_finite() {
                return Err(invalid(
                    "half-plane phase-reference displacement is not representable",
                ));
            }
            let amplitude = component.pressure_pa * cis(phase);
            if !finite(amplitude) {
                return Err(invalid(
                    "half-plane incident amplitude is not representable",
                ));
            }
            components.push((kernel, amplitude));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            frame,
            reference_point_m,
            components,
        })
    }

    pub const fn frame(&self) -> HalfPlaneFrame {
        self.frame
    }

    pub const fn reference_point_m(&self) -> [f64; 3] {
        self.reference_point_m
    }

    pub fn component_count(&self) -> usize {
        self.components.len()
    }

    pub fn acoustic(&self) -> FluidAcousticState {
        self.components[0].0.acoustic()
    }

    /// Original uniform mean-fluid velocity in world coordinates [m/s].
    pub fn mean_flow_velocity_m_s(&self) -> [f64; 3] {
        let speed = self.components[0].0.edge_flow_velocity_m_s;
        self.frame.axes[2].map(|axis| speed * axis)
    }

    /// Coherent total p/u and active/reactive pressure work in the world frame.
    /// Use [`Self::evaluate_energy_transport`] to include mean-flow advection.
    pub fn evaluate(
        &self,
        receiver: HalfPlaneReceiver,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCompactField, AeroacError> {
        let line = self.evaluate_harmonic_line(receiver, cancelled)?;
        let flux = line
            .particle_velocity_m_s
            .map(|u| line.pressure_pa.scale(0.5) * u.conj());
        if flux.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "half-plane angular-spectrum intensity is not representable",
            ));
        }
        Ok(CoherentCompactField {
            pressure_pa: line.pressure_pa,
            particle_velocity_m_s: line.particle_velocity_m_s,
            active_intensity_w_m2: flux.map(|v| v.re),
            reactive_intensity_w_m2: flux.map(|v| v.im),
            scope: if self.components[0].0.edge_flow_velocity_m_s == 0.0 {
                SCOPE
            } else {
                FLOW_SCOPE
            },
        })
    }

    /// Mean perturbation energy and laboratory flux of the complete coherent
    /// spectrum. All cross-component terms survive, including the axial
    /// advection of interference energy, which varies along a flowing edge.
    pub fn evaluate_energy_transport(
        &self,
        receiver: HalfPlaneReceiver,
        cancelled: impl FnMut() -> bool,
    ) -> Result<crate::curle3d::energy::AcousticEnergyTransport, AeroacError> {
        let line = self.evaluate_harmonic_line(receiver, cancelled)?;
        super::energy_transport(self.acoustic(), self.mean_flow_velocity_m_s(), line)
    }

    /// Supply the existing finite-window/sampled observations without a new
    /// delay, conjugation, pressure-to-velocity shortcut or inferred coherence.
    pub fn evaluate_harmonic_line(
        &self,
        receiver: HalfPlaneReceiver,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        poll(&mut cancelled)?;
        let local = self.frame.local_receiver(receiver)?;
        let mut real = [Dd::ZERO; 4];
        let mut imag = [Dd::ZERO; 4];
        for (kernel, amplitude) in &self.components {
            let field = kernel.amplitudes(local, *amplitude, &mut cancelled)?;
            for index in 0..4 {
                real[index] = real[index] + Dd::from_f64(field[index].re);
                imag[index] = imag[index] + Dd::from_f64(field[index].im);
            }
        }
        // Rotate the coherent sum once; retain the DD low words until the
        // complete world component is formed, including near cancellation.
        let velocity = core::array::from_fn(|world| {
            let mut re = Dd::ZERO;
            let mut im = Dd::ZERO;
            for local in 0..3 {
                let basis = Dd::from_f64(self.frame.axes[local][world]);
                re = re + real[local + 1] * basis;
                im = im + imag[local + 1] * basis;
            }
            C64::new(re.to_f64(), im.to_f64())
        });
        let pressure = C64::new(real[0].to_f64(), imag[0].to_f64());
        if !finite(pressure) || velocity.iter().any(|&v| !finite(v)) {
            return Err(invalid("half-plane coherent field is not representable"));
        }
        poll(&mut cancelled)?;
        Ok(HarmonicFieldLine {
            frequency_hz: self.components[0].0.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
        })
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> Dd {
    let mut result = Dd::ZERO;
    for axis in 0..3 {
        result = result + Dd::from_f64(a[axis]) * Dd::from_f64(b[axis]);
    }
    result
}

fn translated_dot(direction: [f64; 3], point: [f64; 3], origin: [f64; 3]) -> Dd {
    let mut result = Dd::ZERO;
    for axis in 0..3 {
        result = result
            + Dd::from_f64(direction[axis])
                * (Dd::from_f64(point[axis]) - Dd::from_f64(origin[axis]));
    }
    result
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|axis| {
        let i = (axis + 1) % 3;
        let j = (axis + 2) % 3;
        (Dd::from_f64(a[i]) * Dd::from_f64(b[j]) - Dd::from_f64(a[j]) * Dd::from_f64(b[i])).to_f64()
    })
}

fn normalized(vector: [f64; 3]) -> Result<[f64; 3], AeroacError> {
    if vector.iter().any(|v| !v.is_finite()) {
        return Err(invalid("half-plane frame directions must be finite"));
    }
    let scale = vector.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    if scale == 0.0 {
        return Err(invalid("half-plane frame directions must be nonzero"));
    }
    let scaled = vector.map(|v| v / scale);
    let length = det::hypot(det::hypot(scaled[0], scaled[1]), scaled[2]);
    Ok(scaled.map(|v| v / length))
}

#[cfg(test)]
mod tests;
