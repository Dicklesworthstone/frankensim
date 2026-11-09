//! One incoming isothermal acoustic wave on the original regularized outlet.
//!
//! Characteristic boundary framework: Poinsot & Lele, JCP 101 (1992),
//! 104–129, <https://doi.org/10.1016/0021-9991(92)90046-2>. The map below is
//! only its linear isothermal plane-normal specialization; the existing
//! regularized lattice reconstruction supplies the microscopic boundary.

use super::*;

/// Linear plane-normal acoustic control about a fixed isothermal reference.
///
/// The pressure face's held density is `rho_ref`; the native lattice supplies
/// `c = sqrt(CS2)` and `p = CS2*rho`. With OUTWARD normal velocity `u_n`,
/// `J+ = (u_n-u_ref) + c*(rho-rho_ref)/rho_ref` leaves the domain and
/// `J- = (u_n-u_ref) - c*(rho-rho_ref)/rho_ref` enters it. This control holds
/// only J- and takes J+ and the tangential velocity from the streamed donor.
/// A zero incoming perturbation admits the reference throughflow, not rest.
///
/// This is the local linear, plane-normal characteristic approximation, not
/// nonlinear NSCBC or a full thermodynamic gas boundary. Reference, donor and
/// reconstructed normal speeds must be strictly outward and subsonic. Backflow
/// needs external tangential/inflow data and refuses; it is never clamped.
/// The caller's density-departure envelope bounds use of the linearization.
/// Native collision/Mach bounds remain independently authoritative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CharacteristicOutlet3 {
    /// Positive outward reference speed [lattice lengths per lattice step].
    pub reference_normal_velocity_lu: f64,
    /// Authored incoming perturbation J- [lattice speed]; zero is undisturbed.
    pub incoming_characteristic_lu: f64,
    /// Inclusive bound on `abs(rho/rho_ref - 1)` at donor and reconstructed
    /// trace. Finite in (0, 1); small values are appropriate to linear acoustics.
    pub maximum_relative_density_departure: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::CollisionModel3;
    use crate::lattice3::{D3q19, GridBudget3};

    #[test]
    fn characteristic_outlet_retains_actual_normal_gradient_for_rheology_and_sgs() {
        for face in Face3::ALL {
            let axis = face as usize / 2;
            let sign = if face as usize % 2 == 0 { -1.0 } else { 1.0 };
            let mut periodic = [true; 3];
            periodic[axis] = false;
            let mut inlet = [0.0; 3];
            inlet[axis] = sign * 0.04;
            let grid = Grid3::<D3q19, 19>::try_new_equilibrium(
                [5; 3],
                CollisionModel3::Bgk { tau: 0.9 },
                [0.0; 3],
                periodic,
                GridBudget3 {
                    maximum_cells: 125,
                    maximum_mach: 0.2,
                },
                |_| false,
                |_| (1.3, inlet),
            )
            .unwrap();
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[face as usize] = FaceBoundary3::Pressure { density: 1.3 };
            faces[face as usize ^ 1] = FaceBoundary3::Velocity { velocity: inlet };
            let mut flow =
                OpenBoundaryGrid3::try_new_with_poll(grid, BoundarySpec3::new(faces), &mut || {
                    Ok(())
                })
                .unwrap();
            let slope = [0.001, -0.002, 0.003];
            let velocity: Vec<_> = (0..125)
                .map(|cell| {
                    let point = coordinates(cell, [5; 3]);
                    core::array::from_fn(|a| inlet[a] + slope[a] * point[axis] as f64)
                })
                .collect();
            let mut point = [2; 3];
            point[axis] = if sign < 0.0 { 0 } else { 4 };
            let cell = (point[2] * 5 + point[1]) * 5 + point[0];
            for row in flow.rheology_gradient(cell, &velocity) {
                assert_eq!(row[axis], 0.0);
            }
            flow.set_characteristic_outlet(
                face,
                CharacteristicOutlet3 {
                    reference_normal_velocity_lu: 0.04,
                    incoming_characteristic_lu: 0.0,
                    maximum_relative_density_departure: 0.05,
                },
                &CancelGate::new_clock_free(),
            )
            .unwrap();
            let gradient = flow.rheology_gradient(cell, &velocity);
            for a in 0..3 {
                assert!((gradient[a][axis] - slope[a]).abs() < 2e-17);
            }
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct PreparedOutlet3 {
    control: CharacteristicOutlet3,
    sound_speed: f64,
}

impl PreparedOutlet3 {
    fn new<E>(control: CharacteristicOutlet3, cs2: f64) -> Result<Self, GridError3<E>> {
        let sound_speed = fs_math::det::sqrt(cs2);
        if !sound_speed.is_finite()
            || sound_speed <= 0.0
            || !control.reference_normal_velocity_lu.is_finite()
            || control.reference_normal_velocity_lu <= 0.0
            || control.reference_normal_velocity_lu >= sound_speed
            || !control.incoming_characteristic_lu.is_finite()
            || !control.maximum_relative_density_departure.is_finite()
            || control.maximum_relative_density_departure <= 0.0
            || control.maximum_relative_density_departure >= 1.0
        {
            return Err(GridError3::Input(
                "invalid linear characteristic outlet control",
            ));
        }
        Ok(Self {
            control,
            sound_speed,
        })
    }

    pub(super) fn trace<E>(
        self,
        face: Face3,
        reference_density: f64,
        donor_density: f64,
        donor_velocity: [f64; 3],
        cell: usize,
    ) -> Result<(f64, [f64; 3]), GridError3<E>> {
        let axis = face as usize / 2;
        let sign = if face as usize % 2 == 0 { -1.0 } else { 1.0 };
        let normal = sign * donor_velocity[axis];
        let departure = (donor_density - reference_density) / reference_density;
        if !departure.is_finite()
            || departure.abs() > self.control.maximum_relative_density_departure
        {
            return Err(GridError3::MaterialField {
                cell,
                reason: "characteristic donor exceeds the relative-density envelope",
            });
        }
        if !normal.is_finite() || normal <= 0.0 || normal >= self.sound_speed {
            return Err(GridError3::MaterialField {
                cell,
                reason: "characteristic outlet requires subsonic outward donor flow",
            });
        }
        let outgoing =
            normal - self.control.reference_normal_velocity_lu + self.sound_speed * departure;
        let incoming = self.control.incoming_characteristic_lu;
        let normal = self.control.reference_normal_velocity_lu + 0.5 * (outgoing + incoming);
        let departure = 0.5 * (outgoing - incoming) / self.sound_speed;
        let density = reference_density + reference_density * departure;
        if !density.is_finite()
            || density <= 0.0
            || !departure.is_finite()
            || departure.abs() > self.control.maximum_relative_density_departure
        {
            return Err(GridError3::MaterialField {
                cell,
                reason: "characteristic trace exceeds the relative-density envelope",
            });
        }
        if !normal.is_finite() || normal <= 0.0 || normal >= self.sound_speed {
            return Err(GridError3::MaterialField {
                cell,
                reason: "characteristic outlet requires subsonic outward reconstructed flow",
            });
        }
        let mut velocity = donor_velocity;
        velocity[axis] = sign * normal;
        Ok((density, velocity))
    }
}

impl<L: EquilibriumLattice3<Q>, const Q: usize> OpenBoundaryGrid3<L, Q> {
    /// Attach or replace a held characteristic control on an existing pressure
    /// face without changing populations, force history or accepted time.
    ///
    /// The face's uniform density becomes the fixed reference. A pressure
    /// profile conflicts with this uniform-reference model; clear it first.
    /// Only complete fluid open nodes are admitted: wall/open seams would
    /// constrain populations independently and cannot impose both acoustic
    /// moments. Interior obstacles are otherwise unchanged. Changing the
    /// reference density requires clearing this control before rebinding.
    ///
    /// Stepping reuses the native donor's physical half-force velocity and
    /// stress, native equilibrium/source transform and measured reservoir
    /// transfer. Serial, pooled, force, rheology/SGS and thermal paths share
    /// that one reconstruction. No callback, history buffer or allocation is
    /// added to the step. No oblique transparency or reflection bound is claimed.
    /// All admission/cancellation failures preserve the previous control.
    pub fn set_characteristic_outlet(
        &mut self,
        face: Face3,
        control: CharacteristicOutlet3,
        gate: &CancelGate,
    ) -> Result<(), GridError3<L::Error>> {
        let poll = || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        let FaceBoundary3::Pressure { density } = self.plan.boundaries.face(face) else {
            return Err(GridError3::Input(
                "characteristic outlet requires a pressure face",
            ));
        };
        if self.pressure_profiles[face as usize].is_some() {
            return Err(GridError3::Input(
                "clear the pressure profile before attaching a characteristic outlet",
            ));
        }
        let prepared = PreparedOutlet3::new(control, L::CS2)?;
        let axis = face as usize / 2;
        let sign = if face as usize % 2 == 0 { -1.0 } else { 1.0 };
        let mut reference = [0.0; 3];
        reference[axis] = sign * control.reference_normal_velocity_lu;
        initial::admit_velocity::<L, Q>(
            self.grid.collision,
            reference,
            0,
            self.grid.budget.maximum_mach,
        )?;
        let (_, target) = prepared.trace(face, density, density, reference, 0)?;
        initial::admit_velocity::<L, Q>(
            self.grid.collision,
            target,
            0,
            self.grid.budget.maximum_mach,
        )?;
        let data = StepData::<L, Q> {
            dims: self.grid.dims,
            periodic: self.grid.periodic,
            solid: &self.grid.solid,
            collision: self.grid.collision,
            collisions: None,
            force: self.grid.force,
            forces: None,
            budget: self.grid.budget,
            marker: PhantomData,
        };
        let mut count = 0;
        for (ordinal, cell) in velocity::plane_cells(self.grid.dims, face).enumerate() {
            if ordinal.is_multiple_of(64) {
                poll()?;
            }
            if self.grid.solid[cell] {
                continue;
            }
            let (walls, _) = self.plan.masks(&data, cell);
            if walls != 0 {
                return Err(GridError3::MaterialField {
                    cell,
                    reason: "characteristic outlet cannot share wall-owned populations",
                });
            }
            count += 1;
        }
        if count == 0 {
            return Err(GridError3::Input(
                "characteristic outlet has no fluid nodes",
            ));
        }
        poll()?;
        self.characteristic_outlets[face as usize] = Some(prepared);
        Ok(())
    }

    /// Restore the original uniform pressure reservoir, without a fluid step.
    pub fn clear_characteristic_outlet(&mut self, face: Face3) {
        self.characteristic_outlets[face as usize] = None;
    }

    /// The accepted incoming-wave control. Its reference density is the held
    /// pressure value in [`Self::boundaries`]. This observation is pure.
    pub fn characteristic_outlet(&self, face: Face3) -> Option<CharacteristicOutlet3> {
        self.characteristic_outlets[face as usize].map(|prepared| prepared.control)
    }
}
