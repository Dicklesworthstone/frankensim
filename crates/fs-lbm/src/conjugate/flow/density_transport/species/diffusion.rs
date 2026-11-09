//! Barycentric binary Fick diffusion on accepted, fixed carrier inventories.
//!
//! One scalar law j_A = -rho D grad(Y_A), j_B = -j_A transports both original
//! constituent masses. D is explicitly prescribed, never inferred from viscosity.
//! Half-cell resistances give one shared flux across a spatial coefficient jump.
//! Each face has an exact positive exponential flow; forward and reverse half
//! sweeps are a symmetric second-order splitting for frozen linear coefficients.
//! A large stable step is not an accuracy certificate: refine physical substeps.
//!
//! This is an isothermal fixed-carrier constitutive component. It does not evolve
//! partial enthalpy, temperature, pressure, chemical potentials or mixture volume.
//! Nonideal thermodynamic factors, Soret/barodiffusion and more than two species
//! require their own original constitutive owners.
//!
//! Fick mass-frame convention: COMSOL Chemical Reaction Engineering, eq. 3-86,
//! https://doc.comsol.com/6.3/doc/com.comsol.help.chem/chem_ug_chemsptrans.08.241.html.
//! Symmetric splitting: Strang (1968), https://doi.org/10.1137/0705041.

use super::{CancelGate, ChtError, MaterialAdvectionStep, Sum, nonnegative, poll, poll_every, ratio};
use crate::conjugate::VoxelDomain;
use fs_math::det;

/// Explicit diffusion work, storage and independent constituent-mass budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinarySpeciesDiffusionControls {
    /// Equal physical substeps, each a forward/reverse pair of half sweeps.
    /// Positive; refine this count to resolve splitting error.
    pub substeps: usize,
    /// Bound on two times the active face count times substeps.
    pub maximum_pair_updates: usize,
    /// Conservative payload bound for cell scratch, output and prepared faces.
    /// Existing carrier storage and allocator overhead are additional.
    pub maximum_storage_bytes: usize,
    /// Per-species local inventory-representation defects plus the absolute
    /// closed-domain mass change [kg]. Trace species retain their own allowance.
    pub maximum_mass_defect_kg: [f64; 2],
    /// Sum of absolute initial and final per-cell partition defects [kg],
    /// measured against the SAME immutable accepted carrier mass.
    pub maximum_partition_defect_kg: f64,
}

/// Physical readback of one complete closed-domain diffusion interval.
#[derive(Clone, Debug, PartialEq)]
pub struct BinarySpeciesDiffusionReport {
    pub duration_s: f64,
    pub substeps: usize,
    /// Fluid-fluid faces with positive molecular conductance, including faces
    /// whose advective mass flux is exactly zero.
    pub diffusive_faces: usize,
    pub pair_updates: usize,
    /// Largest 1-exp(-G*dt_half*(1/M_l+1/M_r)); this measures stiffness,
    /// not the temporal error of a sweep.
    pub maximum_pair_relaxation: f64,
    pub mass_before_kg: [f64; 2],
    pub mass_after_kg: [f64; 2],
    pub mass_change_kg: [f64; 2],
    pub absolute_mass_defect_kg: [f64; 2],
    pub initial_partition_defect_kg: f64,
    pub endpoint_partition_defect_kg: f64,
}

/// Both authoritative constituent inventories, with no carrier or input mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct BinarySpeciesDiffusionStep {
    pub species_mass_kg: Vec<[f64; 2]>,
    pub report: BinarySpeciesDiffusionReport,
}

#[derive(Clone, Copy)]
struct Face {
    cells: [usize; 2],
    /// H*(1-exp(-G*dt_half/H)), H = M_l*M_r/(M_l+M_r) [kg].
    exchange_capacity_kg: f64,
}

#[derive(Clone, Copy, Default)]
struct Inventory {
    mass: [Sum; 2],
    increment: [Sum; 2],
}

fn refuse(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.binary_species_diffusion",
        reason: reason.into(),
    }
}

fn positive(value: f64) -> Result<f64, ChtError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(refuse("positive diffusion geometry or coefficient is unrepresentable"))
    }
}

impl MaterialAdvectionStep {
    /// Diffuse two constituent masses on this step's accepted carrier field.
    ///
    /// Supply post-advection masses, e.g. from advect_species::<2>, or the
    /// previous diffusion endpoint. Their sum must partition this SAME accepted
    /// carrier mass. The original public carrier mass is checked for edits first.
    /// The domain must be the original accepted geometry; its cell count and
    /// fluid mask are checked, while identical physical widths remain the caller's
    /// obligation as in the original pressure/material path.
    ///
    /// Each fluid cell supplies one nonnegative D [m²/s], held fixed throughout
    /// this interval. Density is the accepted mass divided by native cell volume.
    /// For neighboring cells, G=A/(h_l/(2*rho_l*D_l)+h_r/(2*rho_r*D_r)) [kg/s].
    /// This series resistance supports graded Cartesian cells and spatial
    /// coefficient jumps. Zero D blocks the adjacent face; solids and all
    /// exterior faces are impermeable. No boundary concentration is invented.
    ///
    /// Paired transfers update BOTH constituent inventories; the second is never
    /// reconstructed as carrier minus first, and neither is clipped/normalized.
    /// Each pair approaches its mass-weighted equilibrium without overshoot.
    /// Reusable prepared rates keep exponentials outside the sweep loop; work is
    /// O(cells + faces*substeps), storage O(cells + faces), independent of stiffness.
    /// Physical substep and mesh refinement remain necessary for accuracy.
    ///
    /// No thermodynamic heat account is supplied by anonymous species indices.
    /// This does not promote an authored D into a sourced material law or turn
    /// binary Fick transport into arbitrary multicomponent diffusion.
    /// Failure/cancellation publishes no partial state and changes no input.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_binary_species(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; 2]],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        controls: BinarySpeciesDiffusionControls,
        gate: &CancelGate,
    ) -> Result<BinarySpeciesDiffusionStep, ChtError> {
        self.diffuse_binary_species_checked(
            domain, species_mass_kg, diffusivity_m2_s, duration_s, controls, gate,
            |_| poll(gate),
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn diffuse_binary_species_checked(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; 2]],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        controls: BinarySpeciesDiffusionControls,
        gate: &CancelGate,
        mut checkpoint: impl FnMut(usize) -> Result<(), ChtError>,
    ) -> Result<BinarySpeciesDiffusionStep, ChtError> {
        poll(gate)?;
        let transport = self.accepted_mass_transport(gate)?;
        let count = domain.cell_count();
        if count != self.mass_kg.len()
            || species_mass_kg.len() != count
            || diffusivity_m2_s.len() != count
            || controls.substeps == 0
            || controls.substeps as u64 > (1u64 << 53)
        {
            return Err(refuse("binary diffusion dimensions or substeps are invalid"));
        }
        positive(duration_s)?;
        for value in controls.maximum_mass_defect_kg
            .into_iter().chain([controls.maximum_partition_defect_kg])
        {
            nonnegative(value)?;
        }
        let half_dt = positive(0.5 * (duration_s / controls.substeps as f64))?;
        let mut initial_partition = Sum::default();
        let mut face_count = 0usize;
        for c in 0..count {
            poll_every(gate, c)?;
            nonnegative(diffusivity_m2_s[c])?;
            for &mass in &species_mass_kg[c] {
                nonnegative(mass)?;
            }
            if domain.is_fluid(c) != (transport.row_of[c] != usize::MAX) {
                return Err(refuse("diffusion fluid mask differs from the accepted carrier"));
            }
            if !domain.is_fluid(c) {
                if species_mass_kg[c] != [0.0; 2] || diffusivity_m2_s[c] != 0.0 {
                    return Err(refuse("solid cells cannot carry binary species or diffusivity"));
                }
                continue;
            }
            positive(domain.volume(c))?;
            let mut partition = Sum::default();
            partition.add(-self.mass_kg[c]);
            for &mass in &species_mass_kg[c] {
                ratio(mass, self.mass_kg[c])?;
                partition.add(mass);
            }
            initial_partition.add(partition.value.abs());
            for side in [1, 3, 5] {
                if let Some(next) = domain.neighbor(c, side)
                    && domain.is_fluid(next)
                    && diffusivity_m2_s[c] > 0.0
                    && diffusivity_m2_s[next] > 0.0
                {
                    face_count = face_count.checked_add(1)
                        .ok_or_else(|| refuse("diffusion face extent overflows"))?;
                }
            }
        }
        nonnegative(initial_partition.value)?;
        if initial_partition.value > controls.maximum_partition_defect_kg {
            return Err(refuse("binary species do not partition the accepted carrier mass"));
        }
        let pair_updates = face_count.checked_mul(2)
            .and_then(|n| n.checked_mul(controls.substeps))
            .ok_or_else(|| refuse("diffusion work extent overflows"))?;
        if pair_updates > controls.maximum_pair_updates {
            return Err(refuse("binary diffusion pair-update budget exhausted"));
        }
        let bytes = count.checked_mul(core::mem::size_of::<Inventory>() + 16)
            .and_then(|n| n.checked_add(face_count.checked_mul(core::mem::size_of::<Face>())?))
            .and_then(|n| n.checked_add(2048))
            .ok_or_else(|| refuse("diffusion storage extent overflows"))?;
        if bytes > controls.maximum_storage_bytes {
            return Err(refuse("binary diffusion storage budget exhausted"));
        }
        let mut faces = Vec::with_capacity(face_count);
        let mut maximum_pair_relaxation: f64 = 0.0;
        for c in 0..count {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) || diffusivity_m2_s[c] == 0.0 {
                continue;
            }
            for side in [1, 3, 5] {
                let Some(next) = domain.neighbor(c, side).filter(|&n| domain.is_fluid(n))
                else { continue };
                if diffusivity_m2_s[next] == 0.0 {
                    continue;
                }
                let axis = side / 2;
                let rho_left = positive(self.mass_kg[c] / domain.volume(c))?;
                let rho_right = positive(self.mass_kg[next] / domain.volume(next))?;
                let left = positive((0.5 * domain.widths(c)[axis])
                    / positive(rho_left * diffusivity_m2_s[c])?)?;
                let right = positive((0.5 * domain.widths(next)[axis])
                    / positive(rho_right * diffusivity_m2_s[next])?)?;
                let largest = left.max(right);
                let conductance = positive((domain.face_area(c, axis) / largest)
                    / (1.0 + left.min(right) / largest))?;
                let small = self.mass_kg[c].min(self.mass_kg[next]);
                let large = self.mass_kg[c].max(self.mass_kg[next]);
                let capacity = positive(small / (1.0 + small / large))?;
                let mut exponent = (conductance / capacity) * half_dt;
                if !exponent.is_normal() {
                    // An infinite rate-time product means full pair equilibration;
                    // balancing the factors also retains finite scaled products.
                    exponent = det::exp(det::ln(conductance) - det::ln(capacity) + det::ln(half_dt));
                }
                let relaxation = -det::expm1(-exponent);
                if !relaxation.is_finite() || !(0.0..=1.0).contains(&relaxation) {
                    return Err(refuse("binary diffusion relaxation is unrepresentable"));
                }
                maximum_pair_relaxation = maximum_pair_relaxation.max(relaxation);
                faces.push(Face {
                    cells: [c, next],
                    exchange_capacity_kg: positive(capacity * relaxation)?,
                });
            }
        }
        let mut inventories = Vec::with_capacity(count);
        for (c, &mass) in species_mass_kg.iter().enumerate() {
            poll_every(gate, c)?;
            inventories.push(Inventory {
                mass: mass.map(|value| Sum { value, compensation: 0.0 }),
                increment: [Sum::default(); 2],
            });
        }
        // With no conducting faces, no repeated empty sweeps consume a huge
        // caller-specified substep count.
        if !faces.is_empty() {
            for substep in 0..controls.substeps {
                poll(gate)?;
                for (index, face) in faces.iter().enumerate() {
                    poll_every(gate, index)?;
                    exchange(face, &self.mass_kg, &mut inventories)?;
                }
                for (index, face) in faces.iter().rev().enumerate() {
                    poll_every(gate, index)?;
                    exchange(face, &self.mass_kg, &mut inventories)?;
                }
                checkpoint(substep + 1)?;
            }
        }
        let mut before = [Sum::default(); 2];
        let mut after = [Sum::default(); 2];
        let mut change = [Sum::default(); 2];
        let mut defect = [Sum::default(); 2];
        let mut endpoint_partition = Sum::default();
        let mut result = Vec::with_capacity(count);
        for c in 0..count {
            poll_every(gate, c)?;
            let endpoint = inventories[c].mass.map(|m| m.value);
            let mut partition = Sum::default();
            partition.add(-self.mass_kg[c]);
            for k in 0..2 {
                nonnegative(endpoint[k])?;
                before[k].add(species_mass_kg[c][k]);
                after[k].add(endpoint[k]);
                let delta = endpoint[k] - species_mass_kg[c][k];
                change[k].add(delta);
                defect[k].add((delta - inventories[c].increment[k].value).abs());
                partition.add(endpoint[k]);
            }
            endpoint_partition.add(partition.value.abs());
            result.push(endpoint);
        }
        let absolute_mass_defect_kg = [0, 1].map(|k| defect[k].value + change[k].value.abs());
        for k in 0..2 {
            for value in [before[k].value, after[k].value, absolute_mass_defect_kg[k]] {
                nonnegative(value)?;
            }
            if absolute_mass_defect_kg[k] > controls.maximum_mass_defect_kg[k] {
                return Err(refuse("binary diffusion exceeds its constituent-mass budget"));
            }
        }
        let partition_defect = initial_partition.value + endpoint_partition.value;
        nonnegative(partition_defect)?;
        if partition_defect > controls.maximum_partition_defect_kg {
            return Err(refuse("binary diffusion does not preserve the accepted mass partition"));
        }
        poll(gate)?;
        Ok(BinarySpeciesDiffusionStep {
            species_mass_kg: result,
            report: BinarySpeciesDiffusionReport {
                duration_s,
                substeps: controls.substeps,
                diffusive_faces: face_count,
                pair_updates,
                maximum_pair_relaxation,
                mass_before_kg: before.map(|s| s.value),
                mass_after_kg: after.map(|s| s.value),
                mass_change_kg: change.map(|s| s.value),
                absolute_mass_defect_kg,
                initial_partition_defect_kg: initial_partition.value,
                endpoint_partition_defect_kg: endpoint_partition.value,
            },
        })
    }
}

fn exchange(face: &Face, carrier: &[f64], inventory: &mut [Inventory]) -> Result<(), ChtError> {
    let [left, right] = face.cells;
    let a_left = ratio(inventory[left].mass[0].value, carrier[left])?;
    let a_right = ratio(inventory[right].mass[0].value, carrier[right])?;
    let b_left = ratio(inventory[left].mass[1].value, carrier[left])?;
    let b_right = ratio(inventory[right].mass[1].value, carrier[right])?;
    // Complementary gradients are mathematically identical. Use the less
    // abundant constituent so rounded majority fractions cannot erase or
    // overwhelm a trace gradient.
    let gap = if a_left + a_right <= b_left + b_right {
        a_left - a_right
    } else {
        b_right - b_left
    };
    let transferred = face.exchange_capacity_kg * gap;
    if !transferred.is_finite() || (gap != 0.0 && transferred == 0.0) {
        return Err(refuse("nonzero binary diffusion transfer is unrepresentable"));
    }
    for (cell, transfer) in [(left, -transferred), (right, transferred)] {
        for (k, delta) in [transfer, -transfer].into_iter().enumerate() {
            inventory[cell].mass[k].add(delta);
            inventory[cell].increment[k].add(delta);
            nonnegative(inventory[cell].mass[k].value)?;
            if !inventory[cell].mass[k].compensation.is_finite()
                || !inventory[cell].increment[k].value.is_finite()
                || !inventory[cell].increment[k].compensation.is_finite()
            {
                return Err(refuse("binary diffusion inventory arithmetic overflowed"));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
