//! Backward-Euler inertia on the original P1 inelastic Newton solve.
//!
//! The same endpoint constitutive stress, tangent, committed history and
//! backtracking owner solve `M (v_next-v_old)/dt + f_internal = f_external`,
//! with `v_next = (u_next-u_old)/dt`. Reference density supplies consistent
//! mass through the elastic owner's shared P1 element kernel. No material
//! response, geometry, constraints or return mapping is replaced here.
//!
//! Endpoint forces and support reactions do work through the actual increment.
//! Their inertial part equals the kinetic change plus `1/2 delta_v M delta_v`.
//! The latter is time-discretization loss, never physical material heat.
//! Constitutive storage/dissipation and thermal feedback remain with their
//! existing owners. Refine the physical timestep to control this first-order
//! approximation; a work identity is not a temporal-accuracy certificate.

pub mod fatigue;

use super::*;

/// Inertia's contribution to the complete mechanical/thermal work balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InertialEnergy {
    /// Original closed-matter inventory; it is not updated from endpoint T.
    pub reference_mass_kg: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Evaluated as `average(v) dot M delta(v)` to avoid subtracting large K.
    pub kinetic_change_j: f64,
    /// Nonnegative backward-Euler kinetic loss, distinct from plastic loss.
    pub numerical_loss_j: f64,
    pub inertial_work_j: f64,
    /// Inertial work minus kinetic change and kinetic numerical loss.
    pub balance_error_j: f64,
}

/// Complete proposed motion; the original state/history is never mutated.
#[derive(Clone, Debug)]
pub struct BackwardEulerStep<S> {
    /// Native stress, displacement, plastic history and full support reactions.
    pub solution: SmallStrainTetSolution<S>,
    pub velocity_mps: Vec<[f64; 3]>,
    pub energy: InertialEnergy,
}

fn invalid(what: &str) -> SmallStrainTetError {
    SmallStrainTetError::InvalidInput { what: what.into() }
}
fn finite(value: f64) -> Result<f64, SmallStrainTetError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable inertial force, tangent or work"))
    }
}

pub(super) struct Inertia<'a> {
    mass: Csr,
    reference_mass_kg: f64,
    old_u: &'a [f64],
    old_v: &'a [f64],
    inverse_dt: f64,
}
impl<'a> Inertia<'a> {
    fn new(
        elements: &[Element],
        densities: &[f64],
        old_u: &'a [f64],
        old_v: &'a [f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, SmallStrainTetError> {
        if densities.len() != elements.len()
            || densities.iter().any(|rho| !rho.is_finite() || *rho <= 0.0)
            || old_v.len() != old_u.len()
            || old_v.iter().any(|v| !v.is_finite())
            || !dt_s.is_finite()
            || dt_s <= 0.0
        {
            return Err(invalid(
                "one positive reference density per tet, finite full velocity and positive dt required",
            ));
        }
        let inverse_dt = finite(dt_s.recip())?;
        let coefficient = finite(inverse_dt * inverse_dt)?;
        if coefficient <= 0.0 {
            return Err(invalid(
                "inertial tangent coefficient must be representable and positive",
            ));
        }
        let mut mass = Coo::new(old_u.len(), old_u.len());
        let mut reference_mass_kg = 0.0;
        for (element, &density) in elements.iter().zip(densities) {
            checkpoint(cx)?;
            let local = crate::linear3::p1_tet_consistent_mass(element.volume, density);
            let element_mass = finite(density * element.volume)?;
            if element_mass <= 0.0 {
                return Err(invalid(
                    "each tetrahedron needs representable positive reference mass",
                ));
            }
            reference_mass_kg = finite(reference_mass_kg + element_mass)?;
            for a in 0..4 {
                for b in 0..4 {
                    for c in 0..3 {
                        let value = finite(local[3 * a + c][3 * b + c])?;
                        if value <= 0.0 {
                            return Err(invalid("consistent element mass underflowed"));
                        }
                        mass.push(3 * element.nodes[a] + c, 3 * element.nodes[b] + c, value);
                    }
                }
            }
        }
        let mass = mass.assemble();
        for row in 0..old_u.len() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            if !mass.get(row, row).is_finite()
                || mass.get(row, row) <= 0.0
                || mass.row(row).1.iter().any(|m| !m.is_finite())
            {
                return Err(invalid(
                    "every nodal component needs finite positive reference inertia",
                ));
            }
        }
        Ok(Self {
            mass,
            reference_mass_kg,
            old_u,
            old_v,
            inverse_dt,
        })
    }

    // Normalize only physical forces on free rows. Elastic stiffness times
    // rigid travel is unrelated to momentum balance and can accept v_next=0.
    // The positive mass action retains a reference under momentum cancellation;
    // unrelated prescribed-support forces cannot inflate the stopping gate.
    pub(super) fn force_scale(
        &self,
        u: &[f64],
        free: &[usize],
        elastic_free_norm: f64,
        external_free_norm: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, SmallStrainTetError> {
        let mut acceleration_scale = Vec::with_capacity(u.len());
        for (i, ((u, old_u), old_v)) in u.iter().zip(self.old_u).zip(self.old_v).enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            acceleration_scale.push(finite(
                (finite(u - old_u)?.abs() * self.inverse_dt + old_v.abs()) * self.inverse_dt,
            )?);
        }
        // Every entry of the original P1 consistent mass is nonnegative,
        // so its existing action is exactly |M| times this positive vector.
        let force = self.action(&acceleration_scale, cx)?;
        let inertial_norm = free.iter().fold(0.0_f64, |n, &d| n.hypot(force[d]));
        Ok(finite(inertial_norm)?
            .max(finite(elastic_free_norm)?)
            .max(finite(external_free_norm)?)
            .max(f64::MIN_POSITIVE))
    }

    fn action(&self, x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, SmallStrainTetError> {
        let mut output = Vec::with_capacity(x.len());
        for row in 0..x.len() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = self.mass.row(row);
            output.push(finite(
                columns
                    .iter()
                    .zip(values)
                    .fold(0.0, |sum, (&column, &value)| value.mul_add(x[column], sum)),
            )?);
        }
        Ok(output)
    }

    pub(super) fn add_force(
        &self,
        u: &[f64],
        force: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), SmallStrainTetError> {
        let acceleration = u
            .iter()
            .zip(self.old_u)
            .zip(self.old_v)
            .map(|((&u, &old_u), &old_v)| {
                finite(((u - old_u) * self.inverse_dt - old_v) * self.inverse_dt)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let inertial = self.action(&acceleration, cx)?;
        for (force, inertia) in force.iter_mut().zip(inertial) {
            *force = finite(*force + inertia)?;
        }
        Ok(())
    }

    pub(super) fn add_tangent(
        &self,
        matrix: &mut Coo,
        cx: &Cx<'_>,
    ) -> Result<(), SmallStrainTetError> {
        let coefficient = self.inverse_dt * self.inverse_dt;
        for row in 0..self.mass.nrows() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = self.mass.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                let scaled = finite(coefficient * value)?;
                if scaled <= 0.0 {
                    return Err(invalid("scaled consistent inertia underflowed"));
                }
                matrix.push(row, column, scaled);
            }
        }
        Ok(())
    }

    fn motion(
        &self,
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, InertialEnergy), SmallStrainTetError> {
        let du: Vec<_> = u
            .iter()
            .zip(self.old_u)
            .map(|(u, old)| finite(u - old))
            .collect::<Result<_, _>>()?;
        let velocity: Vec<_> = du
            .iter()
            .map(|du| finite(du * self.inverse_dt))
            .collect::<Result<_, _>>()?;
        let delta: Vec<_> = velocity
            .iter()
            .zip(self.old_v)
            .map(|(v, old)| finite(v - old))
            .collect::<Result<_, _>>()?;
        let mean: Vec<_> = velocity
            .iter()
            .zip(self.old_v)
            .map(|(v, old)| finite(0.5 * v + 0.5 * old))
            .collect::<Result<_, _>>()?;
        let old_momentum = self.action(self.old_v, cx)?;
        let next_momentum = self.action(&velocity, cx)?;
        let impulse = self.action(&delta, cx)?;
        let kinetic_before_j = finite(0.5 * dot(self.old_v, &old_momentum))?;
        let kinetic_after_j = finite(0.5 * dot(&velocity, &next_momentum))?;
        let kinetic_change_j = finite(dot(&mean, &impulse))?;
        let numerical_loss_j = finite(0.5 * dot(&delta, &impulse))?;
        let mut inertia = vec![0.0; u.len()];
        self.add_force(u, &mut inertia, cx)?;
        let inertial_work_j = finite(dot(&du, &inertia))?;
        let balance_error_j = finite(inertial_work_j - kinetic_change_j - numerical_loss_j)?;
        let scale = kinetic_before_j
            .max(kinetic_after_j)
            .max(inertial_work_j.abs())
            .max(numerical_loss_j);
        if kinetic_before_j < 0.0
            || kinetic_after_j < 0.0
            || numerical_loss_j < 0.0
            || balance_error_j.abs() > 2048.0 * f64::EPSILON * scale
        {
            return Err(invalid(
                "consistent inertia misses its discrete kinetic work identity",
            ));
        }
        Ok((
            velocity.as_chunks::<3>().0.to_vec(),
            InertialEnergy {
                reference_mass_kg: self.reference_mass_kg,
                kinetic_before_j,
                kinetic_after_j,
                kinetic_change_j,
                numerical_loss_j,
                inertial_work_j,
                balance_error_j,
            },
        ))
    }
}

impl<L: SmallStrainLaw> SmallStrainTetProblem<'_, L> {
    /// One physical inertial interval from the accepted full displacement,
    /// velocity and element history. `load_path` must contain only `1.0`:
    /// forces and prescribed displacements are actual endpoint values.
    /// Support reactions include consistent-mass coupling to constrained DOFs.
    /// Free bodies need no artificial restraints.
    ///
    /// Densities belong to the original reference inventory, in tet order.
    /// A material/thermal caller must retain that inventory across intervals;
    /// changing it requires a separate momentum/matter-transfer law. Every
    /// Newton proposal uses the same old constitutive history. This method
    /// reports inertial work but does not relabel it as constitutive heat.
    /// A motion increment too small to resolve against accumulated absolute
    /// displacement can refuse convergence; force tolerances are not relaxed.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_backward_euler(
        &self,
        displacement_m: &[f64],
        velocity_mps: &[[f64; 3]],
        element_states: &[L::State],
        reference_density_kg_m3: &[f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<BackwardEulerStep<L::State>, SmallStrainTetError> {
        checkpoint(cx)?;
        if self.load_path != [1.0] {
            return Err(invalid(
                "one endpoint load at factor one is required per inertial interval",
            ));
        }
        let elements = self.validate()?;
        self.check(displacement_m, element_states)?;
        let inertia = Inertia::new(
            &elements,
            reference_density_kg_m3,
            displacement_m,
            velocity_mps.as_flattened(),
            dt_s,
            cx,
        )?;
        let solution =
            self.solve_candidate(Some((displacement_m, element_states)), Some(&inertia), cx)?;
        let (velocity_mps, energy) = inertia.motion(solution.displacement_m.as_flattened(), cx)?;
        checkpoint(cx)?;
        Ok(BackwardEulerStep {
            solution,
            velocity_mps,
            energy,
        })
    }
}
