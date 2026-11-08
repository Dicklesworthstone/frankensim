//! Full material-shell motion and the original interior pressure field.
//!
//! The native shell evaluates its complete tensor-Prony/Kelvin operator at
//! one frequency without factoring the bare structure. The native interior
//! owner then factors ALL original displacement, rotation and pressure-potential
//! coordinates together. A bare structural or acoustic pole is not inverted.
//! The fluid volume sets the interface normal; actual P1 translation traces
//! and their exact transpose exchange pressure work once. Physical rotations
//! retain their original inertia and mechanical coupling.
//!
//! This is fixed-geometry, small-signal harmonic coupling in exp(-i omega t).
//! The material is evaluated again for each frequency, and full-cycle strain
//! admission precedes publication. The frequency-bound owner deliberately
//! exposes no time realization or frozen-coefficient frequency derivative.

use super::{
    PreparedStructuralInterior, StructuralInterior, StructuralInteriorError,
    StructuralInteriorResponse, StructuralMatrices, admit, poll,
};
use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::heterogeneous::time::ResistiveBoundary;
use fs_feec::acoustics::helmholtz::heterogeneous::{HeterogeneousHelmholtzFe, VelocityTrace};
use fs_math::{c64::C64, dd::Dd};
use fs_modal::harmonic::HarmonicResponseControls;
use fs_plate::shell::viscosity::harmonic::{
    PreparedTensorPronyShellOperator, TensorPronyShell, TensorShellHarmonicError,
};
use std::collections::BTreeSet;

/// Corresponding original vertices on one planar fluid face and shell facet.
/// The order of either side may be reversed together. Winding never sets the
/// physical normal; the adjacent fluid volume does. No offset skin is inferred.
#[derive(Clone, Copy, Debug)]
pub struct ConformingShellFace {
    pub fluid: [u32; 3],
    pub shell: [usize; 3],
}

/// Source, complete coupled solve and additional interface admission.
#[derive(Clone, Copy, Debug)]
pub struct TensorShellInteriorControls {
    /// Original material assembly and complete shell strain recovery.
    pub material: HarmonicResponseControls,
    /// Counts ALL free shell coordinates plus ALL original pressure coordinates.
    pub coupled: HarmonicResponseControls,
    pub maximum_interface_faces: usize,
    /// Additional topology/trace scratch, conservatively
    /// 128*T + 1024*I + 64*q + 4096 bytes. T is original shell facets,
    /// I selected interfaces, q free shell coordinates. The original geometry,
    /// material operators, coupled pencil and allocator metadata are separate.
    pub maximum_trace_bytes: usize,
    /// Topology/normal-trace construction score only:
    /// 512*(T+1)*(1+ceil(log2(T+1))) + 2048*I + 128*q.
    /// Original material and coupled FE assembly retain their native admission.
    pub maximum_trace_work: u64,
    /// Relative structural, acoustic and enclosing power defects in (0,1).
    pub power_relative_tolerance: f64,
}

#[derive(Debug)]
pub enum TensorShellInteriorError {
    Invalid(&'static str),
    Interior(StructuralInteriorError),
    Shell(TensorShellHarmonicError),
    Budget {
        quantity: &'static str,
        required: u128,
        maximum: u128,
    },
    PowerBalance {
        structure: f64,
        acoustic: f64,
        total: f64,
        tolerance: f64,
    },
}
impl core::fmt::Display for TensorShellInteriorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tensor shell interior: {self:?}")
    }
}
impl std::error::Error for TensorShellInteriorError {}
impl From<StructuralInteriorError> for TensorShellInteriorError {
    fn from(value: StructuralInteriorError) -> Self {
        Self::Interior(value)
    }
}
impl From<TensorShellHarmonicError> for TensorShellInteriorError {
    fn from(value: TensorShellHarmonicError) -> Self {
        Self::Shell(value)
    }
}

/// One original material frequency and its complete coupled spatial pencil.
/// Repeated force/flow cases share a prepared factor. Preparing a different
/// frequency requires the source material owner to evaluate its law again.
pub struct TensorShellInterior<'shell, 'a, 'f, 'c> {
    material: PreparedTensorPronyShellOperator<'shell>,
    model: StructuralInterior<'a, 'f, 'c>,
    tolerance: f64,
}

impl<'shell, 'a, 'f, 'c> TensorShellInterior<'shell, 'a, 'f, 'c> {
    /// Bind selected original shell facets to coincident planar cavity faces.
    /// The full heterogeneous high-order pressure basis, including constant
    /// compression, survives. Unselected fluid faces retain their native rigid
    /// or explicitly supplied resistive boundary. An empty interface is an
    /// explicit detached limit, not an inferred wetted surface.
    ///
    /// Each selected shell facet and fluid face must appear exactly once.
    /// Supported translations have zero trace; rotations are not promoted to
    /// fictitious normal velocities at this zero-offset reference surface.
    /// A curved face requires its own conforming physical trace and refuses.
    /// Resistive admittance is relative to the moving wall, as in the native
    /// interior owner. No additional cavity spring or radiation debit is added.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        interior: &'a HeterogeneousHelmholtzFe<'f, 'c>,
        shell: &'shell TensorPronyShell,
        interfaces: &[ConformingShellFace],
        walls: &[ResistiveBoundary],
        omega_rad_s: f64,
        controls: TensorShellInteriorControls,
        cx: &Cx<'_>,
    ) -> Result<Self, TensorShellInteriorError> {
        poll(cx)?;
        let q = shell.model().free;
        admit(q, interior.ndof(), controls.coupled)?;
        cap(
            "shell material coordinates",
            q as u128,
            controls.material.maximum_coordinates as u128,
        )?;
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(TensorShellInteriorError::Invalid(
                "finite positive material frequency required",
            ));
        }
        if !controls.power_relative_tolerance.is_finite()
            || controls.power_relative_tolerance <= 0.0
            || controls.power_relative_tolerance >= 1.0
        {
            return Err(TensorShellInteriorError::Invalid(
                "finite power tolerance in (0,1) required",
            ));
        }
        let (t, i, n) = (
            shell.mesh().tris.len() as u128,
            interfaces.len() as u128,
            q as u128,
        );
        cap(
            "shell interface faces",
            i,
            controls.maximum_interface_faces as u128,
        )?;
        cap(
            "shell interface trace bytes",
            128 * t + 1024 * i + 64 * n + 4096,
            controls.maximum_trace_bytes as u128,
        )?;
        let levels = if t == 0 {
            0
        } else {
            u128::BITS - t.leading_zeros()
        };
        cap(
            "shell interface trace work",
            512 * (t + 1) * (1 + u128::from(levels)) + 2048 * i + 128 * n,
            controls.maximum_trace_work as u128,
        )?;
        let mut original_facets = BTreeSet::new();
        for &facet in &shell.mesh().tris {
            poll(cx)?;
            let mut key = facet;
            key.sort_unstable();
            if !original_facets.insert(key) {
                return Err(TensorShellInteriorError::Invalid(
                    "shell contains duplicate original facets",
                ));
            }
        }
        let mut selected_shell = BTreeSet::new();
        let mut selected_fluid = BTreeSet::new();
        let mut traces = vec![Vec::new(); q];
        for face in interfaces {
            poll(cx)?;
            let (mut shell_key, mut fluid_key) = (face.shell, face.fluid);
            shell_key.sort_unstable();
            fluid_key.sort_unstable();
            if !original_facets.contains(&shell_key)
                || !selected_shell.insert(shell_key)
                || !selected_fluid.insert(fluid_key)
            {
                return Err(TensorShellInteriorError::Invalid(
                    "interface must pair unique original shell facets and fluid faces",
                ));
            }
            let (points, normal, _) = interior
                .geometry()
                .boundary_geometry(face.fluid)
                .map_err(StructuralInteriorError::from)?;
            for local in 0..3 {
                let node = face.shell[local];
                if points[local] != shell.mesh().nodes[node] {
                    return Err(TensorShellInteriorError::Invalid(
                        "fluid and shell interface vertices must coincide exactly",
                    ));
                }
                for axis in 0..3 {
                    if let Some(coordinate) = shell.model().dof_map[6 * node + axis] {
                        if normal[axis] != 0.0 {
                            let mut values = [0.0; 3];
                            values[local] = normal[axis];
                            traces[coordinate].push((face.fluid, values));
                        }
                    }
                }
            }
        }
        let material = shell.prepare_harmonic_operator(omega_rad_s, controls.material, cx)?;
        let model = StructuralInterior::new(
            interior,
            StructuralMatrices {
                mass: material.mass(),
                stiffness: material.stiffness(),
                damping: Some(material.damping()),
            },
            &traces,
            walls,
            controls.coupled,
            cx,
        )?;
        poll(cx)?;
        Ok(Self {
            material,
            model,
            tolerance: controls.power_relative_tolerance,
        })
    }

    pub fn material_operator(&self) -> &PreparedTensorPronyShellOperator<'shell> {
        &self.material
    }
    pub fn interior(&self) -> &'a HeterogeneousHelmholtzFe<'f, 'c> {
        self.model.interior
    }
    pub fn structural_dofs(&self) -> usize {
        self.model.structural_dofs()
    }
    pub fn pressure_dofs(&self) -> usize {
        self.model.pressure_dofs()
    }
    /// Original row-major pressure x free-shell work pairing.
    pub fn coupling(&self) -> &[f64] {
        self.model.coupling()
    }

    /// Native high-order projection of prescribed inward velocity [m/s].
    /// Pass the resulting original pressure-coordinate flow load to `solve`.
    pub fn project_inward_velocity(
        &self,
        traces: &[VelocityTrace],
        cx: &Cx<'_>,
    ) -> Result<Vec<C64>, TensorShellInteriorError> {
        poll(cx)?;
        let load = self.model.project_inward_velocity(traces)?;
        poll(cx)?;
        Ok(load)
    }

    /// Factor the COMPLETE original system at this exact material frequency.
    /// Neither a bare-shell nor a bare-pressure inverse is a prerequisite.
    pub fn prepare(
        &self,
        cx: &Cx<'_>,
    ) -> Result<PreparedTensorShellInterior<'_, 'shell, 'a, 'f, 'c>, TensorShellInteriorError> {
        poll(cx)?;
        let coupled = self
            .model
            .prepare_angular(self.material.omega_rad_s(), cx)?;
        Ok(PreparedTensorShellInterior {
            owner: self,
            coupled,
        })
    }
}

pub struct PreparedTensorShellInterior<'p, 'shell, 'a, 'f, 'c> {
    owner: &'p TensorShellInterior<'shell, 'a, 'f, 'c>,
    coupled: PreparedStructuralInterior<'p, 'a, 'f, 'c>,
}

/// Original physical fields, source/material heat and signed interface work.
/// The native strain-energy account includes the original drilling stiffness;
/// drilling contributes no material heat. No fitted loss or acoustic mode is used.
#[derive(Clone, Debug)]
pub struct TensorShellInteriorResponse {
    /// The original material frequency; every field uses exp(-i omega t).
    pub angular_frequency_rad_s: f64,
    pub field: StructuralInteriorResponse,
    pub apparent_power_scale_w: f64,
    /// Largest structural, acoustic or enclosing relative power defect.
    pub relative_power_balance: f64,
}

impl PreparedTensorShellInterior<'_, '_, '_, '_, '_> {
    /// Solve complete physical force/torque and original inward volume-flow
    /// loads in peak exp(-i omega t) phasors. Both drives may be coherent.
    /// The same pressure performs source and traction work. An excessive
    /// original material strain or work defect returns no accepted field.
    pub fn solve(
        &self,
        structural_force: &[C64],
        inward_flow_load: &[C64],
        cx: &Cx<'_>,
    ) -> Result<TensorShellInteriorResponse, TensorShellInteriorError> {
        poll(cx)?;
        let field = self.coupled.solve(structural_force, inward_flow_load, cx)?;
        self.owner
            .material
            .admit_displacement(&field.structural_displacement, cx)?;
        let mut scale = Dd::from_f64(0.0);
        for ((force, reaction), velocity) in structural_force
            .iter()
            .zip(&field.fluid_force)
            .zip(&field.structural_velocity)
        {
            poll(cx)?;
            scale = scale + Dd::from_f64(0.5 * (force.abs() + reaction.abs()) * velocity.abs());
        }
        for (flow, pressure) in inward_flow_load.iter().zip(&field.pressure_pa) {
            poll(cx)?;
            scale = scale + Dd::from_f64(0.5 * flow.abs() * pressure.abs());
        }
        scale = scale
            + Dd::from_f64(field.average_structural_dissipation_w.abs())
            + Dd::from_f64(field.average_wall_dissipation_w.abs());
        let scale = scale.to_f64();
        let structural = relative(field.structural_power_balance_residual_w().abs(), scale);
        let acoustic = relative(field.acoustic_power_balance_residual_w().abs(), scale);
        let enclosing = field.average_structural_input_power_w
            + field.average_acoustic_input_power_w
            - field.average_structural_dissipation_w
            - field.average_wall_dissipation_w;
        let total = relative(enclosing.abs(), scale)
            .max(relative(field.power_balance_residual_w.abs(), scale));
        let maximum = structural.max(acoustic).max(total);
        if maximum > self.owner.tolerance {
            return Err(TensorShellInteriorError::PowerBalance {
                structure: structural,
                acoustic,
                total,
                tolerance: self.owner.tolerance,
            });
        }
        poll(cx)?;
        Ok(TensorShellInteriorResponse {
            angular_frequency_rad_s: self.owner.material.omega_rad_s(),
            field,
            apparent_power_scale_w: scale,
            relative_power_balance: maximum,
        })
    }
}

fn cap(
    quantity: &'static str,
    required: u128,
    maximum: u128,
) -> Result<(), TensorShellInteriorError> {
    if required > maximum {
        Err(TensorShellInteriorError::Budget {
            quantity,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn relative(defect: f64, scale: f64) -> f64 {
    if !defect.is_finite() || !scale.is_finite() || scale < 0.0 {
        f64::INFINITY
    } else if scale > 0.0 {
        defect / scale
    } else if defect == 0.0 {
        0.0
    } else {
        f64::INFINITY
    }
}
