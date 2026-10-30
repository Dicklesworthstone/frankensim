//! Charge-conserving conduction in a prescribed moving material snapshot.
//!
//! For in-plane velocity `u` and imposed normal flux density `Bz`, the original
//! scalar or symmetric tensor conductivity acts on `E' = -grad(V) + u x B`.
//! The same prepared conductivity Hodge enforces current continuity and the
//! original voltage/current traces. Unloaded exterior boundaries insulate the
//! TOTAL current, including its motional part.
//!
//! Heat is `E'.J`, force is `J x B`, and signed electrical work is `E.J`:
//! `electrical power = irreversible heat + mechanical power`. These quantities
//! share the same cell current, geometry and material source; a motional EMF
//! cannot become an extra Joule source or an unbalanced terminal current.
//!
//! This is a quasistatic fixed-geometry snapshot with authored cell velocities
//! and imposed normal fields. It does not evolve induction, move the mesh,
//! infer material parameters, add Hall response, or integrate fluid/solid
//! motion. Remesh/re-resolve changed geometry or material through its owner.

use super::{
    ConductionBoundaryResponse, ConductionCellResponse, ConductionControls, ConductionCurrentTrace,
    ConductionError, PlanarConduction, finite, gradients, poll, scalar_potential,
};
use fs_material::conductor::ResolvedConductorState;
use fs_material::conductor::planar::PlanarConductorMaterial;
use fs_solver::SolveReport;

/// Prescribed values on each ORIGINAL triangle, in the mesh-coordinate frame.
/// Positive normal flux density points out of the oriented x-y plane.
#[derive(Clone, Copy, Debug)]
pub struct ConductionMotion<'a> {
    pub velocity_m_s: &'a [[f64; 2]],
    pub normal_flux_density_t: &'a [f64],
}

impl<C: PlanarConductorMaterial> PlanarConduction<C> {
    /// Solve charge continuity with `J = Sigma (-grad(V) + u x B)` using the
    /// existing SPD conductivity operator, including rotated anisotropic data.
    /// Velocity/field changes only rebuild the weak load, not the Hodge matrix.
    ///
    /// Boundary conventions and anchoring are identical to
    /// [`Self::solve_with_current_traces`]. An unconnected conductor needs only
    /// its voltage reference trace: the free potential screens a conservative
    /// motional field without an invented closed-circuit braking current.
    pub fn solve_with_motion(
        &self,
        boundary_voltage_v: &[(usize, f64)],
        current_traces: &[ConductionCurrentTrace],
        motion: ConductionMotion<'_>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingConductionField<'_, C>, ConductionError> {
        poll(cancelled)?;
        if !self.contacts().is_empty() {
            return Err(ConductionError::Invalid(
                "moving electrical contacts require an owned surface work and heat law",
            ));
        }
        self.validate_controls(boundary_voltage_v.len(), controls)?;
        let cells = self.mesh.faces().len();
        if motion.velocity_m_s.len() != cells || motion.normal_flux_density_t.len() != cells {
            return Err(ConductionError::Invalid(
                "one prescribed velocity and normal magnetic field per original triangle required",
            ));
        }
        let prescribed =
            scalar_potential::prescribed_boundary(&self.mesh, boundary_voltage_v, cancelled)?;
        let applied = self.current_loads(current_traces, &prescribed, cancelled)?;
        let mut motional_field = Vec::with_capacity(cells);
        let mut motional_load = vec![0.0; self.mesh.vertices().len()];
        let mut has_motional_field = false;
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let [ux, uy] = motion.velocity_m_s[cell];
            let b = motion.normal_flux_density_t[cell];
            finite(ux)?;
            finite(uy)?;
            finite(b)?;
            let m = [finite(uy * b)?, finite(-ux * b)?];
            motional_field.push(m);
            if m == [0.0; 2] {
                continue;
            }
            has_motional_field = true;
            let gradients = gradients(face.map(|v| self.mesh.vertices()[v as usize]))
                .map_err(|_| ConductionError::Invalid("invalid conductor cell gradient"))?;
            let current = self.conductivity[cell].current_density_a_m2(m);
            let j = [finite(current[0])?, finite(current[1])?];
            let volume = self
                .mesh
                .face_measure(cell)
                .ok_or(ConductionError::Invalid("missing conductor cell volume"))?;
            let mut load = [0.0; 3];
            for local in 1..3 {
                load[local] =
                    finite(volume * gradients[local][0].mul_add(j[0], gradients[local][1] * j[1]))?;
            }
            // Partition of unity makes the cell source sum zero. Use the same
            // reference basis as the gauge-invariant electric-field recovery.
            load[0] = finite(-finite(load[1] + load[2])?)?;
            for local in 0..3 {
                let vertex = face[local] as usize;
                motional_load[vertex] = finite(motional_load[vertex] + load[local])?;
            }
        }
        let motion_field = has_motional_field.then_some(motional_field.as_slice());
        let motion_load = has_motional_field.then_some(motional_load.as_slice());
        let mut rhs: Vec<_> = applied.iter().map(|load| load.unwrap_or(0.0)).collect();
        if let Some(load) = motion_load {
            for (right, &motional) in rhs.iter_mut().zip(load) {
                *right = finite(*right + motional)?;
            }
        }
        let (voltage, report) = scalar_potential::solve(
            &self.stiffness,
            &rhs,
            &prescribed,
            controls.linear,
            cancelled,
        )?;
        let response = self.cell_response_with_motion(&voltage, motion_field, cancelled)?;
        let mut conductor_field_v_m = Vec::with_capacity(cells);
        let mut lorentz_force_density_n_m3 = Vec::with_capacity(cells);
        let mut cell_mechanical_power_w = Vec::with_capacity(cells);
        let mut mechanical_power_w = 0.0;
        let mut electrical_power_w = 0.0;
        let mut work_pairing_residual_w = 0.0;
        for (cell, &m) in motional_field.iter().enumerate() {
            poll(cancelled)?;
            let e = response.electric_field_v_m[cell];
            let j = response.current_density_a_m2[cell];
            let [ux, uy] = motion.velocity_m_s[cell];
            let b = motion.normal_flux_density_t[cell];
            let force = [finite(j[1] * b)?, finite(-j[0] * b)?];
            let volume = self
                .mesh
                .face_measure(cell)
                .ok_or(ConductionError::Invalid("missing conductor cell volume"))?;
            let mechanical = finite(volume * ux.mul_add(force[0], uy * force[1]))?;
            let electrical = finite(volume * e[0].mul_add(j[0], e[1] * j[1]))?;
            let residual = finite(finite(response.cell_power_w[cell] + mechanical)? - electrical)?;
            // An L1 check prevents opposite local work errors from cancelling.
            work_pairing_residual_w = finite(work_pairing_residual_w + residual.abs())?;
            mechanical_power_w = finite(mechanical_power_w + mechanical)?;
            electrical_power_w = finite(electrical_power_w + electrical)?;
            conductor_field_v_m.push([finite(e[0] + m[0])?, finite(e[1] + m[1])?]);
            lorentz_force_density_n_m3.push(force);
            cell_mechanical_power_w.push(mechanical);
        }
        let boundary = self.boundary_response(
            &voltage,
            &prescribed,
            &applied,
            motion_load,
            electrical_power_w,
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            cancelled,
        )?;
        let power_balance_residual_w =
            finite(finite(response.power_w + mechanical_power_w)? - boundary.terminal_power_w)?;
        if work_pairing_residual_w > controls.power_tolerance_w
            || power_balance_residual_w.abs() > controls.power_tolerance_w
        {
            return Err(ConductionError::Balance {
                current_residual_a: boundary.boundary_current_a,
                maximum_free_current_residual_a: boundary.maximum_free_current_residual_a,
                power_residual_w: if work_pairing_residual_w > controls.power_tolerance_w {
                    work_pairing_residual_w
                } else {
                    power_balance_residual_w
                },
            });
        }
        poll(cancelled)?;
        Ok(MovingConductionField {
            problem: self,
            voltage_v: voltage,
            velocity_m_s: motion.velocity_m_s.to_vec(),
            normal_flux_density_t: motion.normal_flux_density_t.to_vec(),
            response,
            conductor_field_v_m,
            lorentz_force_density_n_m3,
            cell_mechanical_power_w,
            mechanical_power_w,
            electrical_power_w,
            boundary,
            nodal_applied_current_a: applied,
            work_pairing_residual_w,
            power_balance_residual_w,
            report,
        })
    }
}

/// Accepted current, heat and force with the exact original source/geometry
/// and prescribed velocity/field snapshots. Mechanical power is positive when
/// the electromagnetic force does work on the moving conductor; terminal power
/// is positive into the conductor and can be negative during generation.
pub struct MovingConductionField<'a, C = ResolvedConductorState> {
    problem: &'a PlanarConduction<C>,
    voltage_v: Vec<f64>,
    velocity_m_s: Vec<[f64; 2]>,
    normal_flux_density_t: Vec<f64>,
    response: ConductionCellResponse,
    conductor_field_v_m: Vec<[f64; 2]>,
    lorentz_force_density_n_m3: Vec<[f64; 2]>,
    cell_mechanical_power_w: Vec<f64>,
    mechanical_power_w: f64,
    electrical_power_w: f64,
    boundary: ConductionBoundaryResponse,
    nodal_applied_current_a: Vec<Option<f64>>,
    work_pairing_residual_w: f64,
    power_balance_residual_w: f64,
    report: Option<SolveReport>,
}

impl<C: PlanarConductorMaterial> MovingConductionField<'_, C> {
    pub fn problem(&self) -> &PlanarConduction<C> {
        self.problem
    }
    pub fn voltage_v(&self) -> &[f64] {
        &self.voltage_v
    }
    pub fn velocity_m_s(&self) -> &[[f64; 2]] {
        &self.velocity_m_s
    }
    pub fn normal_flux_density_t(&self) -> &[f64] {
        &self.normal_flux_density_t
    }
    /// Laboratory electric field `E = -grad(V)` on each original cell [V/m].
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        &self.response.electric_field_v_m
    }
    /// Conductor-frame field `E' = E + u x B` used by the conductivity [V/m].
    pub fn conductor_field_v_m(&self) -> &[[f64; 2]] {
        &self.conductor_field_v_m
    }
    pub fn current_density_a_m2(&self) -> &[[f64; 2]] {
        &self.response.current_density_a_m2
    }
    /// Irreversible `E'.J` integrated over each actual extruded cell [W].
    pub fn cell_joule_power_w(&self) -> &[f64] {
        &self.response.cell_power_w
    }
    pub fn joule_power_w(&self) -> f64 {
        self.response.power_w
    }
    /// Original cell `J x B` [N/m3], without interfacial smoothing.
    pub fn lorentz_force_density_n_m3(&self) -> &[[f64; 2]] {
        &self.lorentz_force_density_n_m3
    }
    pub fn cell_mechanical_power_w(&self) -> &[f64] {
        &self.cell_mechanical_power_w
    }
    pub fn mechanical_power_w(&self) -> f64 {
        self.mechanical_power_w
    }
    /// Signed volume integral of the LABORATORY `E.J` [W].
    pub fn electrical_power_w(&self) -> f64 {
        self.electrical_power_w
    }
    /// Independently recovered boundary voltage-current work, positive in [W].
    pub fn terminal_power_w(&self) -> f64 {
        self.boundary.terminal_power_w
    }
    /// Combined voltage reactions and authored inward current loads [A].
    pub fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        &self.boundary.nodal_boundary_current_a
    }
    pub fn nodal_applied_current_a(&self) -> &[Option<f64>] {
        &self.nodal_applied_current_a
    }
    pub fn boundary_current_a(&self) -> f64 {
        self.boundary.boundary_current_a
    }
    pub fn maximum_free_current_residual_a(&self) -> f64 {
        self.boundary.maximum_free_current_residual_a
    }
    /// Sum of absolute local `heat + mechanical - electrical` defects [W].
    pub fn work_pairing_residual_w(&self) -> f64 {
        self.work_pairing_residual_w
    }
    /// Total `heat + mechanical - terminal` defect [W].
    pub fn power_balance_residual_w(&self) -> f64 {
        self.power_balance_residual_w
    }
    pub fn solve_report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
}
