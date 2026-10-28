//! Quasi-static moving-conductor coupling on the existing electrical face graph.
//!
//! Ohm's law is J = sigma(-grad(phi) + u x B). Each conducting half-cell
//! contributes its length times (u x B) to the face voltage. The transpose of
//! that SAME interpolation scatters I dl x B back to mechanical cells. This
//! closes electrode power = conductor-frame heat + Lorentz mechanical power.
//! B is imposed and stationary during the solve; there is no induction equation,
//! magnetization force, Hall term, moving mesh or inferred material coefficient.

use super::{
    ChtError, ElectricSetup, ElectricSolution, ElectricalContact, Node, Voxel,
    VoxelDomain, finite, finite_positive, poll, poll_every, solve_electric_field_with_contact_field,
};
use fs_exec::CancelGate;

mod isolated;

/// Prescribed SI velocity and magnetic flux density on the electrical cells.
#[derive(Clone, Copy, Debug)]
pub struct ElectricMotionInput<'a> {
    /// One velocity [m/s] per cell. Only bulk conductors contribute motional EMF.
    pub velocity_m_s: &'a [[f64; 3]],
    /// One imposed magnetic flux density [T] per cell, in the velocity frame.
    /// The caller owns solenoidal-field validity and the low-induction regime.
    pub magnetic_flux_density_t: &'a [[f64; 3]],
    /// Cell allocation ceiling; checked before material queries and scratch.
    pub max_cells: usize,
    /// Absolute tolerance [W] for electrical/heat/mechanical and transpose-work
    /// balances. This is a numerical residual limit, not a model-error estimate.
    pub power_tolerance_w: f64,
}

/// A frozen-state current solve with conductor-frame heat and magnetic reaction.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingConductorSolution {
    /// Existing electrical field; its Joule term is I^2 R, never J dot E.
    /// Electrode power now includes mechanical power and need not equal heat.
    pub electric: ElectricSolution,
    /// Force density [N/m3] from the face-current interpolation's transpose.
    /// Zero in insulators and ideal electrode cells. Do not add another J x B.
    pub lorentz_force_density_n_m3: Vec<[f64; 3]>,
    /// Sum of cell force dot the supplied velocity [W]; negative for braking.
    pub mechanical_power_w: f64,
    /// Electrode power minus Joule heat minus mechanical power [W].
    pub power_defect_w: f64,
    /// Mechanical power plus the independently summed I * motional voltage [W].
    pub transfer_defect_w: f64,
}

impl ElectricSetup {
    /// Solve moving material-assigned solids or fluids using the SAME current
    /// matrix, electrode constraints, material admission and heat partition as
    /// [`Self::solve_with_cell_contacts`]. Missing resistivity must be an error;
    /// `None` explicitly declares insulation. Callbacks retain their own receipts.
    ///
    /// Each conducting component requires a voltage reference. Unlike stationary
    /// legacy paths, ungrounded conductors are refused rather than discarded:
    /// they can carry circulating motional currents. Contacts have zero thickness
    /// and therefore no motional voltage or mechanical work of their own. Ideal
    /// electrodes are zero-length electrical constraints, not mechanical bodies.
    ///
    /// The returned force/heat are instantaneous inputs to the existing physics
    /// owners, not a time step or a general electromagnetic-field simulation.
    /// Neither input arrays nor domain occupancy are changed on any outcome.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_cell_motion(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        contacts: &[ElectricalContact],
        motion: ElectricMotionInput<'_>,
        resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
        gate: &CancelGate,
    ) -> Result<MovingConductorSolution, ChtError> {
        self.solve_with_cell_motion_topology(domain, temperature, contacts, &[], motion, resolve, gate)
    }

    /// Moving-conductor counterpart with explicitly insulating electrical faces.
    ///
    /// An opening removes that face from current continuity, motional EMF,
    /// Joule heating and reciprocal magnetic force together. Its two bulk cells
    /// remain intact and retain their complete material laws. This models an
    /// insulating slit/coating in the admitted low-induction regime without a
    /// fictitious huge resistance or deleting mechanical/thermal material.
    ///
    /// Open and resistive face lists must be disjoint and unique in either
    /// orientation. Every conducting component still needs a voltage reference;
    /// an opening which leaves an unreferenced moving body refuses rather than
    /// dropping its possible eddy currents. Electrode-free bodies instead use
    /// [`Self::solve_isolated_cell_motion_topology`]. No dielectric charging,
    /// displacement current, arc law or changing-field induction is introduced.
    /// The result remains a frozen-state solve; inputs and topology are borrowed.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_cell_motion_topology(
        &self,
        domain: &VoxelDomain,
        temperature: &[f64],
        contacts: &[ElectricalContact],
        open_faces: &[[usize; 2]],
        motion: ElectricMotionInput<'_>,
        resolve: impl FnMut(usize, Voxel, f64) -> Result<Option<f64>, ChtError>,
        gate: &CancelGate,
    ) -> Result<MovingConductorSolution, ChtError> {
        poll(gate)?;
        let mut transfer = MotionCoupling::new(domain, motion, gate)?;
        let (nodes, rho) = super::contact::resolve_cells(self, domain, temperature, resolve, gate)?;
        let field = super::contact::ContactField::new(domain, &nodes, contacts, gate)?
            .with_open_faces(domain, open_faces, gate)?;
        let electric = solve_electric_field_with_contact_field(
            domain, self, nodes, &rho, field, Some(&mut transfer), gate,
        )?;
        transfer.finish(domain, electric, gate)
    }
}

fn invalid(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput { field: "electric.motion", reason: reason.into() }
}

pub(super) struct MotionCoupling<'a> {
    input: ElectricMotionInput<'a>,
    ux_b: Vec<[f64; 3]>,
    force_n: Vec<[f64; 3]>,
    motional_power_w: f64,
    isolation: Option<isolated::Isolation>,
}

impl<'a> MotionCoupling<'a> {
    fn new(domain: &VoxelDomain, input: ElectricMotionInput<'a>, gate: &CancelGate) -> Result<Self, ChtError> {
        let n = domain.cell_count();
        if n > input.max_cells || input.velocity_m_s.len() != n
            || input.magnetic_flux_density_t.len() != n
            || !input.power_tolerance_w.is_finite() || input.power_tolerance_w < 0.0
        {
            return Err(invalid("motion fields, cell budget or power tolerance are invalid"));
        }
        let mut ux_b = Vec::with_capacity(n);
        for c in 0..n {
            poll_every(gate, c)?;
            finite_positive("electric.motion.cell_volume_m3", domain.volume(c))?;
            let u = input.velocity_m_s[c];
            let b = input.magnetic_flux_density_t[c];
            if !u.iter().chain(&b).all(|v| v.is_finite()) {
                return Err(invalid(format!("nonfinite velocity or magnetic field at cell {c}")));
            }
            let e = cross(u, b);
            if !e.iter().all(|v| v.is_finite()) {
                return Err(invalid(format!("unrepresentable motional field at cell {c}")));
            }
            ux_b.push(e);
        }
        Ok(Self { input, ux_b, force_n: vec![[0.0; 3]; n], motional_power_w: 0.0, isolation: None })
    }

    pub(super) fn potential_gauges(&self) -> Option<&[bool]> {
        self.isolation.as_ref().map(|state| state.potential_gauges.as_slice())
    }

    // c is always a bulk conductor; n is a conductor or an ideal electrode.
    // Match the existing half-cell resistance geometry, including graded grids.
    pub(super) fn emf(&self, domain: &VoxelDomain, c: usize, n: usize, face: usize, neighbor: Node) -> Result<f64, ChtError> {
        let axis = face / 2;
        let lc = 0.5 * domain.widths(c)[axis];
        let ln = if neighbor == Node::Conductor { 0.5 * domain.widths(n)[axis] } else { 0.0 };
        let value = direction(face) * (lc * self.ux_b[c][axis] + ln * self.ux_b[n][axis]);
        finite("electric.motion.face_voltage_v", value)?;
        Ok(value)
    }

    // Exactly once per physical face, with signed current from c to n.
    pub(super) fn deposit(&mut self, domain: &VoxelDomain, c: usize, n: usize, face: usize, neighbor: Node, current_a: f64) -> Result<(), ChtError> {
        let axis = face / 2;
        let mut tangent = [0.0; 3];
        tangent[axis] = direction(face);
        finite("electric.motion.face_current_a", current_a)?;
        if let Some(isolation) = &mut self.isolation {
            isolation.record(c, n, current_a)?;
        }
        self.motional_power_w += current_a * self.emf(domain, c, n, face, neighbor)?;
        for (cell, bulk) in [(c, true), (n, neighbor == Node::Conductor)] {
            if bulk {
                let reaction = cross(tangent, self.input.magnetic_flux_density_t[cell]);
                let weight = current_a * (0.5 * domain.widths(cell)[axis]);
                for (force, value) in self.force_n[cell].iter_mut().zip(reaction) {
                    *force += weight * value;
                    finite("electric.motion.force_n", *force)?;
                }
            }
        }
        finite("electric.motion.motional_power_w", self.motional_power_w)
    }

    fn finish(self, domain: &VoxelDomain, electric: ElectricSolution, gate: &CancelGate) -> Result<MovingConductorSolution, ChtError> {
        if let Some(isolation) = &self.isolation {
            isolation.check(gate)?;
        }
        let mut force_density = Vec::with_capacity(self.force_n.len());
        let mut mechanical_power_w = 0.0;
        for (c, force) in self.force_n.iter().enumerate() {
            poll_every(gate, c)?;
            let mut density = [0.0; 3];
            for a in 0..3 {
                // Work is measured from the actual returned force density,
                // independently of the face-motional-voltage reduction.
                density[a] = force[a] / domain.volume(c);
                finite("electric.motion.force_density_n_m3", density[a])?;
                mechanical_power_w += density[a] * domain.volume(c) * self.input.velocity_m_s[c][a];
            }
            force_density.push(density);
        }
        let power_defect_w = electric.delivered_w - electric.total_joule_w - mechanical_power_w;
        let transfer_defect_w = mechanical_power_w + self.motional_power_w;
        for value in [power_defect_w, transfer_defect_w] {
            if !value.is_finite() || value.abs() > self.input.power_tolerance_w {
                return Err(invalid(format!("electrical/mechanical power defect {value} W exceeds tolerance")));
            }
        }
        poll(gate)?;
        Ok(MovingConductorSolution {
            electric, lorentz_force_density_n_m3: force_density,
            mechanical_power_w, power_defect_w, transfer_defect_w,
        })
    }
}

fn direction(face: usize) -> f64 { if face % 2 == 0 { -1.0 } else { 1.0 } }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
}
