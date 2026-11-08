//! Reciprocal plane-wave fluid ports for the original spatial piezoelectric mesh.
//!
//! The uniform guide mode has volume flow U=b^T v, with b the EXACT P1
//! integral of outward normal displacement over selected exterior faces.
//! For a semi-infinite guide of area A, Z=rho*c/A and p=Z*U+2*p_in.
//! The traction is -b*p, so its loading -i*omega*Z*b*b^T enters the same
//! electroelastic solve. Both incident sensing and outgoing actuation retain
//! fluid backreaction; no separately tuned microphone/actuator gain is used.
//!
//! This is an explicit single uniform-mode guide model. The caller declares
//! its maximum wavenumber and amplitude domain. Higher transverse/evanescent
//! modes, guide discontinuities, viscothermal losses and free-space diffraction
//! require their own acoustic operators; a plane port is not a free-field BEM.
//! Source owners retain their own phase, state and source authority.
//! See COMSOL's Port and Impedance theory:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.028.html>
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.023.html>

use super::super::TetPiezoelectricError;
use super::{C64, Result, TetPiezoelectricModel, current, finite, invalid, poll, positive, power};
use fs_exec::Cx;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::det;
use std::collections::{BTreeMap, BTreeSet};

/// Declared single-mode approximation and cold interface-work limits.
#[derive(Clone, Copy, Debug)]
pub struct PlaneWavePortControls {
    /// Bound on retained nonzero mechanical surface coordinates. The harmonic
    /// rank-one stamp has at most the square of this many entries.
    pub maximum_surface_dofs: usize,
    /// Relative plane/normal tolerance, strictly between zero and one.
    pub coplanarity_tolerance: f64,
    /// Maximum admitted omega/c [rad/m], supplied from the actual guide model.
    /// For example, stay below the first transverse cutoff for a rigid guide.
    pub maximum_wavenumber_rad_m: f64,
    /// Maximum peak |p|/(rho*c^2) and |U|/(A*c), strictly below one.
    /// A small value declares the linear acoustic operating envelope.
    pub maximum_relative_amplitude: f64,
}

/// Prepared geometric pressure/flow map, borrowed from one original model.
/// The full guide section may include a rigid part: its area must be at least
/// the selected moving-face area. All moving faces must share one outward plane.
pub struct PlaneWaveAcousticPort<'a> {
    pub(in crate::piezoelectric::spatial) model: &'a TetPiezoelectricModel<'a>,
    pub(in crate::piezoelectric::spatial) faces: Vec<[usize; 3]>,
    weights_m2: Vec<(usize, f64)>,
    moving_area_m2: f64,
    section_area_m2: f64,
    pub(in crate::piezoelectric::spatial) controls: PlaneWavePortControls,
}

impl PlaneWaveAcousticPort<'_> {
    /// Actual moving-face area from the original mesh [m^2].
    pub const fn moving_area_m2(&self) -> f64 {
        self.moving_area_m2
    }
    /// Full acoustic guide section, including any rigid part [m^2].
    pub const fn section_area_m2(&self) -> f64 {
        self.section_area_m2
    }
    /// Original mechanical DOF and integrated outward normal shape weight.
    pub fn weights_m2(&self) -> &[(usize, f64)] {
        &self.weights_m2
    }
}

/// One original fluid owner and incoming peak pressure at its prepared port.
/// Use zero incoming pressure for an anechoic outgoing load. Material state is
/// evaluated once per harmonic solve, before assembly, and never inside Krylov.
#[derive(Clone, Copy)]
pub struct HarmonicAcousticPort<'a> {
    /// Prepared map on the exact model passed to the harmonic solve.
    pub port: &'a PlaneWaveAcousticPort<'a>,
    /// Original immutable fluid/material/EOS owner; no air default is supplied.
    pub medium: &'a dyn FluidAcoustics,
    /// Incoming traveling-wave pressure at this section [Pa], exp(-i omega t).
    pub incident_pressure_pa: C64,
}

/// Solved acoustic waves and power in the same retained uniform mode.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicAcousticPortResponse {
    /// Coherent density/c/bulk-modulus snapshot from the original owner.
    pub medium: FluidAcousticState,
    /// Full guide cross-section [m^2].
    pub section_area_m2: f64,
    /// Volume velocity out of the solid and into the guide [m^3/s].
    pub volume_velocity_m3_s: C64,
    /// Total compressive fluid pressure acting on the solid [Pa].
    pub pressure_pa: C64,
    /// Prescribed incoming traveling-wave pressure [Pa].
    pub incident_pressure_pa: C64,
    /// Computed outgoing traveling-wave pressure [Pa].
    pub outgoing_pressure_pa: C64,
    /// Incoming wave power, nonnegative [W].
    pub average_incident_power_w: f64,
    /// Outgoing wave power, nonnegative [W].
    pub average_outgoing_power_w: f64,
    /// Signed mechanical power entering the guide [W]. Evaluated from p,U
    /// directly to avoid subtracting nearly equal incident/outgoing powers.
    pub average_net_outward_power_w: f64,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

impl TetPiezoelectricModel<'_> {
    /// Reuse an admitted geometric trace when only the material snapshot changes.
    /// Exact borrowed geometry and frame identity make the original exterior
    /// topology, normals and shape integrals authoritative. This avoids another
    /// whole-mesh face-owner scan for every finite-temperature Newton trial.
    pub(in crate::piezoelectric::spatial) fn rebind_acoustic_port(
        &self,
        original: &PlaneWaveAcousticPort<'_>,
        cx: &Cx<'_>,
    ) -> Result<PlaneWaveAcousticPort<'_>> {
        poll(cx)?;
        if !std::ptr::eq(self.nodes_m, original.model.nodes_m)
            || !std::ptr::eq(self.tetrahedra, original.model.tetrahedra)
            || self.mesh_frame != original.model.mesh_frame
        {
            return Err(invalid(
                "acoustic trace refresh requires the exact original geometry and frame",
            ));
        }
        let mut faces = Vec::with_capacity(original.faces.len());
        for face in &original.faces {
            poll(cx)?;
            faces.push(*face);
        }
        let mut weights_m2 = Vec::with_capacity(original.weights_m2.len());
        for weight in &original.weights_m2 {
            poll(cx)?;
            weights_m2.push(*weight);
        }
        Ok(PlaneWaveAcousticPort {
            model: self,
            faces,
            weights_m2,
            moving_area_m2: original.moving_area_m2,
            section_area_m2: original.section_area_m2,
            controls: original.controls,
        })
    }

    /// Prepare one complete exterior-face pressure trace with its adjoint flow
    /// map. Face order/orientation supplied by the caller does not change the
    /// result: the tetrahedron's opposite vertex defines the outward normal.
    /// Repeated, interior, noncoplanar or oppositely directed faces refuse.
    pub fn plane_wave_acoustic_port(
        &self,
        faces: &[[usize; 3]],
        section_area_m2: f64,
        controls: PlaneWavePortControls,
        cx: &Cx<'_>,
    ) -> Result<PlaneWaveAcousticPort<'_>> {
        poll(cx)?;
        if faces.is_empty()
            || faces.len() > self.tetrahedra.len().saturating_mul(4)
            || !positive(section_area_m2)
            || controls.maximum_surface_dofs == 0
            || !positive(controls.coplanarity_tolerance)
            || controls.coplanarity_tolerance >= 1.0
            || !positive(controls.maximum_wavenumber_rad_m)
            || !positive(controls.maximum_relative_amplitude)
            || controls.maximum_relative_amplitude >= 1.0
        {
            return Err(invalid(
                "bounded nonempty planar acoustic trace and finite positive guide controls required",
            ));
        }
        let mut selected = BTreeSet::new();
        for face in faces {
            poll(cx)?;
            let mut face = *face;
            face.sort_unstable();
            if face[2] >= self.nodes_m.len()
                || face[0] == face[1]
                || face[1] == face[2]
                || !selected.insert(face)
            {
                return Err(invalid("acoustic faces must be unique original triangles"));
            }
        }
        let mut owners = BTreeMap::<[usize; 3], (usize, usize)>::new();
        for tet in self.tetrahedra {
            poll(cx)?;
            for opposite in 0..4 {
                let mut face = [0; 3];
                let mut at = 0;
                for (i, &node) in tet.iter().enumerate() {
                    if i != opposite {
                        face[at] = node;
                        at += 1;
                    }
                }
                face.sort_unstable();
                if selected.contains(&face) {
                    let owner = owners.entry(face).or_insert((0, tet[opposite]));
                    owner.0 += 1;
                }
            }
        }
        let first = *selected.first().expect("nonempty selected trace");
        let origin = self.nodes_m[first[0]];
        let mut normal: Option<[f64; 3]> = None;
        let mut weights = BTreeMap::<usize, f64>::new();
        let mut area = 0.0;
        let mut extent = 0.0_f64;
        for face in &selected {
            for &node in face {
                for distance in sub(self.nodes_m[node], origin) {
                    extent = extent.max(distance.abs());
                }
            }
        }
        for face in &selected {
            poll(cx)?;
            let &(count, opposite) = owners
                .get(face)
                .ok_or_else(|| invalid("acoustic trace is absent from the original mesh"))?;
            if count != 1 {
                return Err(invalid("acoustic port requires exterior faces"));
            }
            let [a, b, c] = face.map(|i| self.nodes_m[i]);
            let mut twice_area_normal = cross(sub(b, a), sub(c, a));
            if dot(twice_area_normal, sub(self.nodes_m[opposite], a)) > 0.0 {
                twice_area_normal = twice_area_normal.map(|v| -v);
            }
            let twice_area = det::sqrt(dot(twice_area_normal, twice_area_normal));
            if !positive(twice_area) {
                return Err(invalid("acoustic face area is unrepresentable"));
            }
            let unit = twice_area_normal.map(|v| v / twice_area);
            let reference = *normal.get_or_insert(unit);
            if sub(unit, reference)
                .iter()
                .any(|v| v.abs() > controls.coplanarity_tolerance)
                || [a, b, c].iter().any(|&p| {
                    dot(sub(p, origin), reference).abs() > controls.coplanarity_tolerance * extent
                })
            {
                return Err(invalid(
                    "one acoustic guide requires coplanar consistently outward faces",
                ));
            }
            area += 0.5 * twice_area;
            for &node in face {
                for axis in 0..3 {
                    let weight = twice_area_normal[axis] / 6.0;
                    if weight != 0.0 {
                        *weights.entry(3 * node + axis).or_insert(0.0) += weight;
                    }
                }
            }
        }
        if !positive(area)
            || area > section_area_m2 * (1.0 + controls.coplanarity_tolerance)
            || weights.len() > controls.maximum_surface_dofs
            || weights.values().any(|w| !w.is_finite())
        {
            return Err(invalid(
                "acoustic area or surface-coordinate budget is inconsistent",
            ));
        }
        poll(cx)?;
        Ok(PlaneWaveAcousticPort {
            model: self,
            faces: selected.into_iter().collect(),
            weights_m2: weights.into_iter().collect(),
            moving_area_m2: area,
            section_area_m2,
            controls,
        })
    }
}

pub(super) struct PreparedAcousticPort<'a> {
    port: &'a PlaneWaveAcousticPort<'a>,
    medium: FluidAcousticState,
    impedance: f64,
    incident: C64,
}

pub(super) fn prepare<'a>(
    model: &TetPiezoelectricModel<'_>,
    inputs: &[HarmonicAcousticPort<'a>],
    omega: f64,
    cx: &Cx<'_>,
) -> Result<Vec<PreparedAcousticPort<'a>>> {
    if inputs.len() > model.tetrahedra.len().saturating_mul(4) {
        return Err(invalid(
            "acoustic port count exceeds original exterior-face capacity",
        ));
    }
    let mut faces = BTreeSet::new();
    let mut prepared = Vec::with_capacity(inputs.len());
    for input in inputs {
        poll(cx)?;
        if !core::ptr::eq(model, input.port.model) || !finite(input.incident_pressure_pa) {
            return Err(invalid(
                "acoustic port must retain its original model and finite incident wave",
            ));
        }
        for &face in &input.port.faces {
            if !faces.insert(face) {
                return Err(invalid("acoustic ports cannot load one face twice"));
            }
        }
        let medium = input
            .medium
            .acoustic_state()
            .map_err(TetPiezoelectricError::AcousticMaterial)?;
        let impedance = medium.characteristic_impedance_pa_s_m() / input.port.section_area_m2;
        if !positive(impedance)
            || omega / medium.sound_speed_m_s() > input.port.controls.maximum_wavenumber_rad_m
            || input.incident_pressure_pa.abs() / medium.adiabatic_bulk_modulus_pa()
                > input.port.controls.maximum_relative_amplitude
        {
            return Err(invalid(
                "fluid acoustic port leaves its admitted frequency or amplitude domain",
            ));
        }
        prepared.push(PreparedAcousticPort {
            port: input.port,
            medium,
            impedance,
            incident: input.incident_pressure_pa,
        });
    }
    Ok(prepared)
}

pub(super) fn add_sources(
    ports: &[PreparedAcousticPort<'_>],
    forces: &mut [[C64; 3]],
    cx: &Cx<'_>,
) -> Result<()> {
    for port in ports {
        for &(dof, b) in &port.port.weights_m2 {
            poll(cx)?;
            let f = &mut forces[dof / 3][dof % 3];
            *f = *f - port.incident.scale(2.0 * b);
            if !finite(*f) {
                return Err(invalid("incident acoustic force overflows"));
            }
        }
    }
    Ok(())
}

pub(super) fn stamp(
    ports: &[PreparedAcousticPort<'_>],
    omega: f64,
    add: &mut impl FnMut(usize, usize, C64, f64) -> Result<()>,
    cx: &Cx<'_>,
) -> Result<()> {
    for port in ports {
        for (i, &(row, a)) in port.port.weights_m2.iter().enumerate() {
            poll(cx)?;
            for &(col, b) in &port.port.weights_m2[i..] {
                let value = omega * (port.impedance * a) * b;
                let entry = C64::new(0.0, -value);
                add(row, col, entry, value.abs())?;
                if row != col {
                    add(col, row, entry, value.abs())?;
                }
            }
        }
    }
    Ok(())
}

pub(super) fn recover(
    ports: &[PreparedAcousticPort<'_>],
    omega: f64,
    full: &[C64],
    effective_force: &mut [[C64; 3]],
    power_tolerance_w: f64,
    cx: &Cx<'_>,
) -> Result<Vec<HarmonicAcousticPortResponse>> {
    let mut responses = Vec::with_capacity(ports.len());
    for port in ports {
        poll(cx)?;
        let mut flow = C64::ZERO;
        for &(dof, b) in &port.port.weights_m2 {
            flow = flow + current(full[dof], omega).scale(b);
        }
        let outgoing = port.incident + flow.scale(port.impedance);
        let pressure = outgoing + port.incident;
        let incident_power = 0.5 * port.incident.norm_sq() / port.impedance;
        let outgoing_power = 0.5 * outgoing.norm_sq() / port.impedance;
        let net = power(pressure, flow);
        // Difference of squared wave amplitudes in factored form avoids
        // subtracting two nearly equal large positive wave powers.
        let wave_net = power(outgoing - port.incident, outgoing + port.incident) / port.impedance;
        let amplitude = port.port.controls.maximum_relative_amplitude;
        if [flow, pressure, outgoing].iter().any(|&v| !finite(v))
            || [incident_power, outgoing_power, net, wave_net]
                .iter()
                .any(|v| !v.is_finite())
            || pressure.abs() / port.medium.adiabatic_bulk_modulus_pa() > amplitude
            || outgoing.abs() / port.medium.adiabatic_bulk_modulus_pa() > amplitude
            || flow.abs() / (port.port.section_area_m2 * port.medium.sound_speed_m_s()) > amplitude
        {
            return Err(invalid(
                "solved acoustic waves leave their finite small-signal domain",
            ));
        }
        if (net - wave_net).abs() > power_tolerance_w {
            return Err(TetPiezoelectricError::Balance(
                "fluid traction work differs from outgoing minus incident wave power",
            ));
        }
        for &(dof, b) in &port.port.weights_m2 {
            poll(cx)?;
            let f = &mut effective_force[dof / 3][dof % 3];
            *f = *f - pressure.scale(b);
        }
        responses.push(HarmonicAcousticPortResponse {
            medium: port.medium,
            section_area_m2: port.port.section_area_m2,
            volume_velocity_m3_s: flow,
            pressure_pa: pressure,
            incident_pressure_pa: port.incident,
            outgoing_pressure_pa: outgoing,
            average_incident_power_w: incident_power,
            average_outgoing_power_w: outgoing_power,
            average_net_outward_power_w: net,
        });
    }
    Ok(responses)
}
