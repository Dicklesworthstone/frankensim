//! Prescribed three-dimensional motion in the original steady-current field.
//!
//! `J = Sigma(E + u cross B)`, with `E = -grad(V)`. The motional weak load
//! changes charge continuity; the original conductivity operator stays SPD.
//! Conductor-frame Joule heat and Lorentz output power remain distinct:
//! `integral(E.J) = integral((E + u cross B).J) + integral(u.(J cross B))`.
//!
//! Velocity and magnetic flux density are constant per original tetrahedron.
//! This fixed-domain, nonrelativistic DC approximation prescribes B and motion;
//! it does not evolve induction, finite moving interfaces or mechanical state.

use super::*;

/// Laboratory-frame prescribed fields on the original tetrahedra. A fixed
/// Eulerian mesh requires motion to preserve the represented material/source
/// occupancy, or a separately justified instantaneous quasistatic snapshot.
#[derive(Clone, Copy, Debug)]
pub struct SpatialConductionMotion<'a> {
    pub velocity_m_s: &'a [[f64; 3]],
    pub flux_density_t: &'a [[f64; 3]],
}

pub(super) struct PreparedMotion<'a> {
    pub(super) prescribed: SpatialConductionMotion<'a>,
    pub(super) field_v_m: Vec<[f64; 3]>,
    pub(super) load_a: Vec<f64>,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> Result<[f64; 3]> {
    Ok([
        finite(a[1].mul_add(b[2], -a[2] * b[1]))?,
        finite(a[2].mul_add(b[0], -a[0] * b[2]))?,
        finite(a[0].mul_add(b[1], -a[1] * b[0]))?,
    ])
}

impl<C: SpatialConductorMaterial> SpatialConduction<C> {
    /// Solve the original terminal/current-trace problem with prescribed
    /// `u cross B`. Open components retain their native zero-current gauge:
    /// a conservative motional field is screened by voltage, not an invented
    /// short circuit. Scalar and full rotated conductivity tensors share this
    /// same solve and source ownership.
    ///
    /// Independent charge, local work, Hodge and boundary-power gates precede
    /// publication. Finite-resistance contacts refuse until their own moving
    /// surface work/heat law is supplied; shared-vertex ideal continuity remains
    /// part of the original volume field. The caller owns magnetic consistency,
    /// motion and thermal evolution outside this quasistatic snapshot.
    pub fn solve_with_motion(
        &self,
        boundary: &SpatialConductionBoundary,
        motion: SpatialConductionMotion<'_>,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialMovingConductionField<'_, C>> {
        check(cancelled)?;
        if !self.contacts.is_empty() {
            return Err(SpatialConductionError::Invalid(
                "moving electrical contacts require an owned surface work and heat law",
            ));
        }
        let cells = self.conductors.len();
        if motion.velocity_m_s.len() != cells || motion.flux_density_t.len() != cells {
            return Err(SpatialConductionError::Invalid(
                "one prescribed velocity and magnetic field per original tetrahedron required",
            ));
        }
        let traces = self.prepare_boundary(boundary, controls, cancelled)?;
        let g = &self.geometry;
        let mut prepared = PreparedMotion {
            prescribed: motion,
            field_v_m: Vec::with_capacity(cells),
            load_a: vec![0.0; g.mesh.vertex_count],
        };
        for (cell, tet) in g.mesh.tets.iter().enumerate() {
            check(cancelled)?;
            for value in motion.velocity_m_s[cell]
                .into_iter()
                .chain(motion.flux_density_t[cell])
            {
                finite(value)?;
            }
            let motional = cross(motion.velocity_m_s[cell], motion.flux_density_t[cell])?;
            prepared.field_v_m.push(motional);
            if motional == [0.0; 3] {
                continue;
            }
            let current = self.conductivity[cell].current_density_a_m2(motional);
            for value in current {
                finite(value)?;
            }
            let volume = g.geometry.vol_signed[cell].abs();
            let mut local_load = [0.0; 4];
            for (local, load) in local_load.iter_mut().enumerate().skip(1) {
                *load = finite(volume * dot(g.geometry.grads[cell][local], current))?;
            }
            // Exact partition of unity in the same reference-node convention
            // as electric-field recovery: motion cannot inject net charge.
            local_load[0] = finite(-finite(
                finite(local_load[1] + local_load[2])? + local_load[3],
            )?)?;
            for (&vertex, load) in tet.iter().zip(local_load) {
                let vertex = vertex as usize;
                prepared.load_a[vertex] = finite(prepared.load_a[vertex] + load)?;
            }
        }
        let mut rhs = traces.flux_load_wb.clone();
        for (value, &motional) in rhs.iter_mut().zip(&prepared.load_a) {
            check(cancelled)?;
            if motional != 0.0 {
                *value = finite(*value + motional)?;
            }
        }
        let (relative, report) = solve_scalar_system(
            &self.stiffness,
            &rhs,
            &traces.prescribed,
            controls.linear,
            cancelled,
        )?;
        let (field, power) = self.recover_field_with_motion(
            relative,
            traces,
            Some(&prepared),
            controls.current_tolerance_a,
            controls.power_tolerance_w,
            report,
            cancelled,
        )?;
        let power = power.ok_or(SpatialConductionError::Invalid(
            "missing original motional power recovery",
        ))?;
        check(cancelled)?;
        Ok(SpatialMovingConductionField {
            problem: self,
            velocity_m_s: motion.velocity_m_s.to_vec(),
            flux_density_t: motion.flux_density_t.to_vec(),
            field,
            power,
        })
    }
}

pub(super) struct MotionPower {
    conductor_field_v_m: Vec<[f64; 3]>,
    lorentz_force_density_n_m3: Vec<[f64; 3]>,
    cell_mechanical_power_w: Vec<f64>,
    pub(super) mechanical_power_w: f64,
    pub(super) electrical_power_w: f64,
    pub(super) work_pairing_residual_w: f64,
}

impl MotionPower {
    pub(super) fn new(cells: usize) -> Self {
        Self {
            conductor_field_v_m: Vec::with_capacity(cells),
            lorentz_force_density_n_m3: Vec::with_capacity(cells),
            cell_mechanical_power_w: Vec::with_capacity(cells),
            mechanical_power_w: 0.0,
            electrical_power_w: 0.0,
            work_pairing_residual_w: 0.0,
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn record(
        &mut self,
        motion: SpatialConductionMotion<'_>,
        cell: usize,
        electric: [f64; 3],
        conductor: [f64; 3],
        current: [f64; 3],
        heat: f64,
        volume: f64,
    ) -> Result<()> {
        let force = cross(current, motion.flux_density_t[cell])?;
        let mechanical = finite(volume * dot(motion.velocity_m_s[cell], force))?;
        let electrical = finite(volume * dot(electric, current))?;
        let residual = finite(finite(heat + mechanical)? - electrical)?;
        self.mechanical_power_w = finite(self.mechanical_power_w + mechanical)?;
        self.electrical_power_w = finite(self.electrical_power_w + electrical)?;
        // Opposite cellwise defects cannot hide in one cancelling global sum.
        self.work_pairing_residual_w = finite(self.work_pairing_residual_w + residual.abs())?;
        self.conductor_field_v_m.push(conductor);
        self.lorentz_force_density_n_m3.push(force);
        self.cell_mechanical_power_w.push(mechanical);
        Ok(())
    }
}

/// Accepted original-cell currents and physical power transfers. Positive
/// mechanical power is delivered TO prescribed motion. Terminal power is
/// positive INTO the conductor and reverses sign during loaded generation.
/// Only Joule power belongs in a thermal receiver's irreversible heat source.
pub struct SpatialMovingConductionField<'a, C = ResolvedConductorState> {
    problem: &'a SpatialConduction<C>,
    velocity_m_s: Vec<[f64; 3]>,
    flux_density_t: Vec<[f64; 3]>,
    field: SpatialConductionField,
    power: MotionPower,
}

impl<C: SpatialConductorMaterial> SpatialMovingConductionField<'_, C> {
    pub fn problem(&self) -> &SpatialConduction<C> {
        self.problem
    }
    pub fn voltage_v(&self) -> &[f64] {
        &self.field.voltage_v
    }
    pub fn velocity_m_s(&self) -> &[[f64; 3]] {
        &self.velocity_m_s
    }
    pub fn flux_density_t(&self) -> &[[f64; 3]] {
        &self.flux_density_t
    }
    /// Laboratory `E=-grad(V)` [V/m].
    pub fn electric_field_v_m(&self) -> &[[f64; 3]] {
        &self.field.electric_field_v_m
    }
    /// Conductor-frame `E'=E+u cross B` [V/m].
    pub fn conductor_field_v_m(&self) -> &[[f64; 3]] {
        &self.power.conductor_field_v_m
    }
    pub fn current_density_a_m2(&self) -> &[[f64; 3]] {
        &self.field.current_density_a_m2
    }
    /// Actual tetrahedron-volume integral of `E'.J` [W].
    pub fn cell_joule_power_w(&self) -> &[f64] {
        &self.field.cell_power_w
    }
    pub fn joule_power_w(&self) -> f64 {
        self.field.power_w
    }
    pub fn lorentz_force_density_n_m3(&self) -> &[[f64; 3]] {
        &self.power.lorentz_force_density_n_m3
    }
    pub fn cell_mechanical_power_w(&self) -> &[f64] {
        &self.power.cell_mechanical_power_w
    }
    pub fn mechanical_power_w(&self) -> f64 {
        self.power.mechanical_power_w
    }
    /// Signed volume integral of laboratory `E.J` [W].
    pub fn electrical_power_w(&self) -> f64 {
        self.power.electrical_power_w
    }
    /// Independently recovered boundary voltage-current work [W].
    pub fn terminal_power_w(&self) -> f64 {
        self.field.supplied_power_w
    }
    pub fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        &self.field.nodal_boundary_current_a
    }
    pub fn nodal_applied_current_a(&self) -> &[f64] {
        &self.field.nodal_applied_current_a
    }
    pub fn boundary_current_a(&self) -> f64 {
        self.field.current_residual_a
    }
    pub fn maximum_free_current_residual_a(&self) -> f64 {
        self.field.maximum_free_current_residual_a
    }
    pub fn work_pairing_residual_w(&self) -> f64 {
        self.power.work_pairing_residual_w
    }
    /// Total `Joule + mechanical - terminal` defect [W].
    pub fn power_balance_residual_w(&self) -> f64 {
        self.field.power_residual_w
    }
    /// Laboratory field power minus original weak nodal-current work [W].
    pub fn hodge_power_residual_w(&self) -> f64 {
        self.field.hodge_power_residual_w
    }
    pub fn gauges(&self) -> &[usize] {
        &self.field.gauges
    }
    pub fn solve_report(&self) -> Option<&SolveReport> {
        self.field.linear.as_ref()
    }
}
