//! Reference areal inventory integrated with the native P1 velocity field.
//!
//! A triangle retains its own three corner values, so material-interface jumps
//! are not averaged at shared nodes. Density is per REFERENCE area; stretching
//! or thermal changes of the current geometry cannot create inventory. These
//! areal values do not assert a through-thickness mass profile or infer an
//! elastic law from density. The original resolved sheets keep that ownership.

use super::{MembraneError, MembranePatch, finite_vectors, poll};

/// Two discretizations of the SAME retained reference areal density.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MembraneMassImage {
    /// Positive row sums. Preserves total mass and first moments, while
    /// approximating the kinetic energy of a nonuniform velocity field.
    Lumped,
    /// Exact P1 density times both P1 velocity shapes: degree-three moments.
    /// Retains inter-node inertia and continuum kinetic moments on each cell.
    Consistent,
}

// c is a positive normalized moment, bounded away from zero by 1/60. If
// rho*area overflows while the individual moment remains finite, apply c
// before restoring rho. No raw sum of large corner densities is formed.
fn weighted_mass(area: f64, rho: f64, c: f64) -> Result<f64, MembraneError> {
    let product = rho * area;
    let value = if product.is_finite() {
        product * c
    } else {
        rho * (area * c)
    };
    if !value.is_finite() || value <= 0.0 {
        return Err(MembraneError::Invalid(
            "unrepresentable reference membrane mass moment",
        ));
    }
    Ok(value)
}

/// Positive row sums of the exact P1 reference-density mass on a triangle.
///
/// `area` is its positive reference area [m²], independent of facet winding.
/// `density` contains reference areal density [kg/m²] at its three corners,
/// in the same order as the velocity shapes. Shared vertices may carry different
/// values on adjacent facets; no material-interface averaging occurs here.
/// These are the original membrane load/lumping weights, not nodal samples
/// multiplied by one third of the area. Their sum and first moment integrate
/// the complete linear density. No current geometry or constitutive law enters.
///
/// # Errors
/// Nonpositive/nonfinite area or corner density, or an unrepresentable moment.
pub fn p1_reference_row_masses(area: f64, density: [f64; 3]) -> Result<[f64; 3], MembraneError> {
    if !area.is_finite() || area <= 0.0 {
        return Err(MembraneError::Invalid(
            "positive finite reference triangle area required",
        ));
    }
    if density.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err(MembraneError::Invalid(
            "positive finite corner reference areal density required",
        ));
    }
    if density[0] == density[1] && density[0] == density[2] {
        // Preserve the original resolved-sheet arithmetic exactly.
        let mass = density[0] * area;
        if !mass.is_finite() || mass / 3.0 <= 0.0 {
            return Err(MembraneError::Invalid(
                "unrepresentable reference membrane mass",
            ));
        }
        return Ok([mass / 3.0; 3]);
    }
    let scale = density.into_iter().fold(0.0_f64, f64::max);
    let r = density.map(|x| x / scale);
    let sum = r[0] + r[1] + r[2];
    let mut masses = [0.0; 3];
    for a in 0..3 {
        masses[a] = weighted_mass(area, scale, (sum + r[a]) / 12.0)?;
    }
    Ok(masses)
}

/// Exact scalar P1 translational inertia for one reference triangle [kg].
///
/// The same area/corner convention as [`p1_reference_row_masses`] applies.
/// `M_ab = integral rho0 N_a N_b dA` retains all cubic density/velocity
/// moments. Reuse this same 3-by-3 block independently for each Cartesian
/// component. Supports restrict both matrix axes after this physical block is
/// formed; fixed nodes never donate mass to free nodes. Density is reference
/// areal inventory, so this function never multiplies thickness again.
///
/// # Errors
/// Nonpositive/nonfinite area or corner density, or an unrepresentable moment.
pub fn p1_reference_consistent_mass(
    area: f64,
    density: [f64; 3],
) -> Result<[[f64; 3]; 3], MembraneError> {
    if !area.is_finite() || area <= 0.0 || density.iter().any(|rho| !rho.is_finite() || *rho <= 0.0)
    {
        return Err(MembraneError::Invalid(
            "positive finite reference triangle area and corner areal density required",
        ));
    }
    let scale = density.into_iter().fold(0.0_f64, f64::max);
    let r = density.map(|x| x / scale);
    let sum = r[0] + r[1] + r[2];
    let mut mass = [[0.0; 3]; 3];
    for a in 0..3 {
        mass[a][a] = weighted_mass(area, scale, (sum + 2.0 * r[a]) / 30.0)?;
        for b in a + 1..3 {
            let value = weighted_mass(area, scale, (sum + r[a] + r[b]) / 60.0)?;
            mass[a][b] = value;
            mass[b][a] = value;
        }
    }
    Ok(mass)
}

impl MembranePatch {
    /// Bind positive P1 reference areal density [kg/m²] at each ORIGINAL
    /// triangle's three corners, in its stored connectivity order.
    ///
    /// This explicit inventory replaces the sheets' scalar areal density for
    /// inertia and acceleration loads. Original source objects, elastic laws,
    /// thickness and grain are unchanged; the field does not acquire measured
    /// material authority merely by being admitted here. A shared node may have
    /// different incident values to preserve a material-interface jump.
    ///
    /// The topology budget was admitted at construction. Density coverage,
    /// representability and cancellation are checked before publishing a new
    /// patch. Subsequent deformation never resamples or changes these values.
    pub fn with_reference_areal_density(
        mut self,
        density_kg_m2: Vec<[f64; 3]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MembraneError> {
        poll(&mut cancelled)?;
        if density_kg_m2.len() != self.elements.len() {
            return Err(MembraneError::Invalid(
                "one reference areal-density triple per original triangle required",
            ));
        }
        let mut masses = vec![0.0; self.mesh.nodes.len()];
        for (element, &density) in self.elements.iter().zip(&density_kg_m2) {
            poll(&mut cancelled)?;
            let local = p1_reference_row_masses(element.area, density)?;
            for a in 0..3 {
                masses[element.nodes[a]] += local[a];
            }
        }
        if masses.iter().any(|m| !m.is_finite() || *m <= 0.0)
            || !masses.iter().sum::<f64>().is_finite()
        {
            return Err(MembraneError::Invalid(
                "unrepresentable assembled reference membrane mass",
            ));
        }
        poll(&mut cancelled)?;
        self.nodal_mass_kg = masses;
        self.reference_areal_density_kg_m2 = density_kg_m2;
        Ok(self)
    }

    /// Immutable corner density, in original triangle/connectivity order.
    /// Without an explicit field each triple repeats its original sheet value.
    pub fn reference_corner_areal_density_kg_m2(&self) -> &[[f64; 3]] {
        &self.reference_areal_density_kg_m2
    }

    /// Current areal density at the material corners of each original triangle.
    /// P1 surface deformation has one constant area ratio J per triangle, so
    /// the exact inventory relation is `rho_area_current = rho_area_reference/J`.
    /// Current area is an observation; it never replaces reference inertia.
    /// The original finite-strain geometry/material domain must still admit the
    /// supplied state. No thickness or bulk-density assumption is introduced.
    pub fn current_corner_areal_density_kg_m2(
        &self,
        displacement_m: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        poll(&mut cancelled)?;
        if !finite_vectors(displacement_m, self.mesh.nodes.len()) {
            return Err(MembraneError::Invalid(
                "one finite membrane displacement per node required",
            ));
        }
        let mut density = Vec::with_capacity(self.elements.len());
        for (index, element) in self.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let kinematics = self.kinematics(element, displacement_m);
            self.response(index, kinematics)?;
            let [a, b] = kinematics.gradient;
            let cross = [
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            ];
            let ratio = cross[0].hypot(cross[1]).hypot(cross[2]);
            let current = self.reference_areal_density_kg_m2[index].map(|rho| rho / ratio);
            if current.iter().any(|rho| !rho.is_finite() || *rho <= 0.0) {
                return Err(MembraneError::Invalid(
                    "unrepresentable current membrane areal density",
                ));
            }
            density.push(current);
        }
        poll(&mut cancelled)?;
        Ok(density)
    }

    /// Apply the selected reference inertia independently to all three spatial
    /// components. Velocity produces nodal momentum [kg m/s]; acceleration
    /// produces force [N]. Consistent inertia integrates both density and the
    /// P1 input field together, including their cross terms.
    ///
    /// This observation/load entry allocates its returned O(nodes) vector.
    /// Failure/cancellation publishes no partial output. Time consumers prepare
    /// their own inverse and reuse their integration work arrays.
    pub fn reference_mass_apply(
        &self,
        image: MembraneMassImage,
        values: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, MembraneError> {
        poll(&mut cancelled)?;
        if !finite_vectors(values, self.mesh.nodes.len()) {
            return Err(MembraneError::Invalid(
                "one finite reference inertia vector per node required",
            ));
        }
        let mut output = vec![[0.0; 3]; values.len()];
        match image {
            MembraneMassImage::Lumped => {
                for (node, (&mass, value)) in self.nodal_mass_kg.iter().zip(values).enumerate() {
                    poll(&mut cancelled)?;
                    output[node] = value.map(|v| mass * v);
                }
            }
            MembraneMassImage::Consistent => {
                for (index, element) in self.elements.iter().enumerate() {
                    poll(&mut cancelled)?;
                    let mass = self.element_consistent_mass(index)?;
                    for a in 0..3 {
                        for b in 0..3 {
                            for axis in 0..3 {
                                output[element.nodes[a]][axis] +=
                                    mass[a][b] * values[element.nodes[b]][axis];
                            }
                        }
                    }
                }
            }
        }
        if !finite_vectors(&output, values.len()) {
            return Err(MembraneError::Invalid(
                "unrepresentable reference membrane inertia action",
            ));
        }
        poll(&mut cancelled)?;
        Ok(output)
    }

    pub(super) fn element_consistent_mass(
        &self,
        index: usize,
    ) -> Result<[[f64; 3]; 3], MembraneError> {
        p1_reference_consistent_mass(
            self.elements[index].area,
            self.reference_areal_density_kg_m2[index],
        )
    }

    // Callers own the dense work cap. Both native modes and prepared dynamics
    // use the same element kernel; unsupported coordinates are eliminated only
    // after constructing the original reference inertia contribution.
    pub(super) fn mass_on_dofs(
        &self,
        image: MembraneMassImage,
        dofs: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, MembraneError> {
        poll(cancelled)?;
        let n = dofs.len();
        let entries = n.checked_mul(n).ok_or(MembraneError::Budget(
            "reference membrane mass matrix extent",
        ))?;
        if dofs.iter().any(|&dof| dof >= 3 * self.mesh.nodes.len())
            || dofs.windows(2).any(|pair| pair[0] >= pair[1])
        {
            return Err(MembraneError::Invalid(
                "ordered unique in-range membrane inertia coordinates required",
            ));
        }
        let mut matrix = vec![0.0; entries];
        match image {
            MembraneMassImage::Lumped => {
                for (row, &dof) in dofs.iter().enumerate() {
                    matrix[row * n + row] = self.nodal_mass_kg[dof / 3];
                }
            }
            MembraneMassImage::Consistent => {
                let mut index = vec![usize::MAX; 3 * self.mesh.nodes.len()];
                for (row, &dof) in dofs.iter().enumerate() {
                    index[dof] = row;
                }
                for (element_index, element) in self.elements.iter().enumerate() {
                    poll(cancelled)?;
                    let local = self.element_consistent_mass(element_index)?;
                    for axis in 0..3 {
                        for a in 0..3 {
                            let row = index[3 * element.nodes[a] + axis];
                            if row == usize::MAX {
                                continue;
                            }
                            for b in 0..3 {
                                let column = index[3 * element.nodes[b] + axis];
                                if column != usize::MAX {
                                    matrix[row * n + column] += local[a][b];
                                }
                            }
                        }
                    }
                }
            }
        }
        if matrix.iter().any(|x| !x.is_finite()) {
            return Err(MembraneError::Invalid(
                "unrepresentable reference membrane mass matrix",
            ));
        }
        poll(cancelled)?;
        Ok(matrix)
    }
}
