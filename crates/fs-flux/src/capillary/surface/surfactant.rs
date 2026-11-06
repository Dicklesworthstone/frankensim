//! Isothermal insoluble surfactant on the existing material-triangle surface.
//!
//! A face carries conserved amount N [mol], so material advection is exactly
//! Gamma=N/A. Its Helmholtz free energy is E=A*psi(Gamma); mechanical tension
//! is dE/dA at fixed N, namely sigma=psi-Gamma*psi'. Using sigma*A as storage
//! would erase the mixing free energy and incorrectly treat dilution as an
//! external coefficient supply. The ideal dilute and finite-packing Langmuir
//! laws below retain this distinction. See Barrett, Garcke and Nurnberg (2015),
//! <https://doi.org/10.1051/m2an/2014039>, for the thermodynamic relation and EOS.
//!
//! The finite-deformation force composes the ORIGINAL triangle area discrete
//! gradient with the exact area secant of E. Its analytic endpoint tangent
//! includes the changing concentration. Geometry, fixed topology, source
//! applicability and represented interface count remain the caller's concern.
//! The caller supplies RT in J/mol at a held temperature; no solubility,
//! adsorption, micelles, surface phase transition or thermal heat is inferred.
//! Negative tensions and occupied packing limits refuse rather than clip.
//! Kernels allocate nothing and use deterministic elementary logarithms.

use super::{
    CapillaryPatchResponse, CapillarySurfaceError, Sum, SurfaceEnergyChange,
    TriangleCapillaryResponse, Vector, cross, dot, geometry, sub, triangle_surface_tension,
    triangle_surface_tension_step, triangle_surface_tension_step_tangent,
};
use fs_math::det;

/// Conservative, dissipative exchange on an explicitly authored surface graph.
pub mod diffusion;

/// Explicit isothermal, nonionic insoluble-monolayer equation of state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SurfactantLaw {
    /// Ideal dilute monolayer: sigma=sigma_clean-RT*Gamma. The reference
    /// concentration fixes the chemical/free-energy datum, not the force.
    IdealDilute {
        /// Tension of the clean, explicitly identified interface [N/m].
        clean_tension_n_m: f64,
        /// Product of the molar gas constant and held temperature [J/mol].
        thermal_energy_j_mol: f64,
        /// Positive concentration defining the chemical datum [mol/m²].
        reference_concentration_mol_m2: f64,
    },
    /// Ideal lattice monolayer: sigma=sigma_clean+RT*Gamma_inf*ln(1-theta),
    /// theta=Gamma/Gamma_inf. This represents finite packing, not adsorption
    /// kinetics or a material-specific critical micelle concentration.
    Langmuir {
        /// Tension of the clean, explicitly identified interface [N/m].
        clean_tension_n_m: f64,
        /// Product of the molar gas constant and held temperature [J/mol].
        thermal_energy_j_mol: f64,
        /// Positive maximum monolayer concentration [mol/m²].
        saturation_concentration_mol_m2: f64,
    },
}

/// Fixed material-face inventory; amounts do not follow changing face area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfactantFace {
    /// Explicit interface equation of state and its held physical parameters.
    pub law: SurfactantLaw,
    /// Conserved, nonnegative surface amount [mol].
    pub amount_mol: f64,
}

/// Borrowed face inventory with a common or explicitly heterogeneous law.
#[derive(Clone, Copy, Debug)]
pub enum SurfactantField<'a> {
    /// One law with an independent conserved amount on each material face.
    Uniform {
        law: SurfactantLaw,
        amounts_mol: &'a [f64],
    },
    /// One law/inventory pair in the supplied face-connectivity order.
    PerFace(&'a [SurfactantFace]),
}

/// Thermodynamic and mechanical state at one actual face area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfactantState {
    /// Conserved amount divided by current area [mol/m²].
    pub concentration_mol_m2: f64,
    /// Total isothermal surface Helmholtz free energy [J], with the declared
    /// chemical datum. This is not sigma*A and need not be nonnegative.
    pub free_energy_j: f64,
    /// Nonnegative mechanical surface tension [N/m].
    pub tension_n_m: f64,
    /// Gibbs dilatational elasticity -Gamma*d(sigma)/d(Gamma) [N/m].
    pub dilatational_elasticity_n_m: f64,
    /// Molar chemical potential [J/mol]. `None` is the -infinity limiting
    /// potential at exactly zero amount; the energy there remains finite.
    pub chemical_potential_j_mol: Option<f64>,
}

/// Invalid input, unavailable EOS state or unrepresentable arithmetic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfactantError {
    Surface(CapillarySurfaceError),
    InvalidLaw,
    InvalidAmount,
    InvalidArea,
    PackingLimit,
    NegativeTension,
    Unrepresentable,
}

impl core::fmt::Display for SurfactantError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "insoluble surface surfactant: {self:?}")
    }
}
impl std::error::Error for SurfactantError {}
impl From<CapillarySurfaceError> for SurfactantError {
    fn from(value: CapillarySurfaceError) -> Self {
        Self::Surface(value)
    }
}

impl SurfactantLaw {
    fn parameters(self) -> (f64, f64, f64) {
        match self {
            Self::IdealDilute {
                clean_tension_n_m,
                thermal_energy_j_mol,
                reference_concentration_mol_m2,
            } => (
                clean_tension_n_m,
                thermal_energy_j_mol,
                reference_concentration_mol_m2,
            ),
            Self::Langmuir {
                clean_tension_n_m,
                thermal_energy_j_mol,
                saturation_concentration_mol_m2,
            } => (
                clean_tension_n_m,
                thermal_energy_j_mol,
                saturation_concentration_mol_m2,
            ),
        }
    }

    pub(super) fn validate(self) -> Result<(), SurfactantError> {
        let (clean, thermal, concentration) = self.parameters();
        if [clean, thermal, concentration]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            Err(SurfactantError::InvalidLaw)
        } else {
            Ok(())
        }
    }
}

impl SurfactantField<'_> {
    pub(super) fn validate(self, faces: usize) -> Result<(), SurfactantError> {
        match self {
            Self::Uniform { law, amounts_mol } if amounts_mol.len() == faces => {
                law.validate()?;
                amounts_mol.iter().try_for_each(|&n| amount(n))
            }
            Self::PerFace(values) if values.len() == faces => values.iter().try_for_each(|face| {
                face.law.validate()?;
                amount(face.amount_mol)
            }),
            _ => Err(CapillarySurfaceError::LengthMismatch.into()),
        }
    }

    pub(super) fn at(self, face: usize) -> SurfactantFace {
        match self {
            Self::Uniform { law, amounts_mol } => SurfactantFace {
                law,
                amount_mol: amounts_mol[face],
            },
            Self::PerFace(values) => values[face],
        }
    }
}

impl SurfactantFace {
    /// Evaluate the same free-energy law used by forces, without an area
    /// reference or inferred amount. Zero inventory recovers the clean surface.
    pub fn at_area(self, area_m2: f64) -> Result<SurfactantState, SurfactantError> {
        self.law.validate()?;
        amount(self.amount_mol)?;
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(SurfactantError::InvalidArea);
        }
        let (clean, rt, reference) = self.law.parameters();
        let clean_energy = clean * area_m2;
        if !clean_energy.is_finite() {
            return Err(SurfactantError::Unrepresentable);
        }
        if self.amount_mol == 0.0 {
            return Ok(SurfactantState {
                concentration_mol_m2: 0.0,
                free_energy_j: clean_energy,
                tension_n_m: clean,
                dilatational_elasticity_n_m: 0.0,
                chemical_potential_j_mol: None,
            });
        }
        let concentration = self.amount_mol / area_m2;
        if !concentration.is_finite() || concentration <= 0.0 {
            return Err(SurfactantError::Unrepresentable);
        }
        let osmotic = rt * concentration;
        let logarithm = log_ratio(concentration, reference);
        let (tension, elasticity, chemical, mixing) = match self.law {
            SurfactantLaw::IdealDilute { .. } => {
                (clean - osmotic, osmotic, rt * logarithm, logarithm - 1.0)
            }
            SurfactantLaw::Langmuir { .. } => {
                let theta = concentration / reference;
                if !theta.is_finite() || theta >= 1.0 {
                    return Err(SurfactantError::PackingLimit);
                }
                let log_vacancy = log1p(-theta);
                (
                    clean - osmotic * log1p_ratio(-theta),
                    osmotic / (1.0 - theta),
                    rt * (logarithm - log_vacancy),
                    logarithm - 1.0 - h(-theta),
                )
            }
        };
        if tension < 0.0 {
            return Err(SurfactantError::NegativeTension);
        }
        let energy = (rt * self.amount_mol).mul_add(mixing, clean_energy);
        if [tension, elasticity, chemical, energy]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(SurfactantError::Unrepresentable);
        }
        Ok(SurfactantState {
            concentration_mol_m2: concentration,
            free_energy_j: energy,
            tension_n_m: tension,
            dilatational_elasticity_n_m: elasticity,
            chemical_potential_j_mol: Some(chemical),
        })
    }

    /// Free-energy change for amount transport on a held face area [J].
    /// The clean area energy cancels algebraically. Both endpoint states must
    /// remain inside the same EOS; no reservoir work or heat is inferred.
    pub fn energy_change_at_fixed_area(
        self,
        area_m2: f64,
        new_amount_mol: f64,
    ) -> Result<f64, SurfactantError> {
        self.at_area(area_m2)?;
        let other = Self {
            amount_mol: new_amount_mol,
            ..self
        };
        other.at_area(area_m2)?;
        if self.amount_mol == new_amount_mol {
            return Ok(0.0);
        }
        let (high, low, sign) = if self.amount_mol > new_amount_mol {
            (self.amount_mol, new_amount_mol, 1.0)
        } else {
            (new_amount_mol, self.amount_mol, -1.0)
        };
        let (_, rt, reference) = self.law.parameters();
        let concentration = high / area_m2;
        let logarithm = log_ratio(concentration, reference);
        let dn = low - high;
        let ratio = dn / high;
        let coefficient = match self.law {
            SurfactantLaw::IdealDilute { .. } => logarithm + h(ratio),
            SurfactantLaw::Langmuir { .. } => {
                let theta = concentration / reference;
                let vacant_ratio = -(dn / area_m2 / reference) / (1.0 - theta);
                logarithm - log1p(-theta) + h(ratio) - h(vacant_ratio)
            }
        };
        let change = sign * (rt * dn) * coefficient;
        if !change.is_finite() {
            return Err(SurfactantError::Unrepresentable);
        }
        Ok(change)
    }
}

/// Conservative finite-deformation force at fixed material-face amount.
/// The area secant supplies BOTH normal capillarity and assembled tangential
/// Marangoni stress from spatially varying inventory. No coefficient supply is
/// added: concentration changes are internal to the same stored free energy.
pub fn triangle_surfactant_step(
    before_m: [Vector; 3],
    after_m: [Vector; 3],
    face: SurfactantFace,
) -> Result<TriangleCapillaryResponse, SurfactantError> {
    let before = geometry(before_m)?;
    let after = geometry(after_m)?;
    let initial = face.at_area(before.area)?;
    let final_state = face.at_area(after.area)?;
    if face.amount_mol == 0.0 {
        return Ok(triangle_surface_tension_step(
            before_m,
            after_m,
            initial.tension_n_m,
            initial.tension_n_m,
        )?);
    }
    let (secant, _) = area_secant(face, before.area, after.area, initial, final_state)?;
    let mut response = triangle_surface_tension_step(before_m, after_m, secant, secant)?;
    response.energy.energy_before_j = initial.free_energy_j;
    response.energy.energy_after_j = final_state.free_energy_j;
    response.energy.energy_change_j = secant * (after.area - before.area);
    response.energy.coefficient_supply_j = 0.0;
    if !response.energy.energy_change_j.is_finite() {
        return Err(SurfactantError::Unrepresentable);
    }
    Ok(response)
}

/// Exact final-position force derivative at held initial geometry and amount.
/// This composes the original area tangent with the derivative of the EOS
/// secant; at equal endpoint areas its scalar derivative is E_AA/2, not zero.
pub fn triangle_surfactant_step_tangent(
    before_m: [Vector; 3],
    after_m: [Vector; 3],
    face: SurfactantFace,
    direction_m: [Vector; 3],
) -> Result<[Vector; 3], SurfactantError> {
    if direction_m.iter().flatten().any(|v| !v.is_finite()) {
        return Err(CapillarySurfaceError::NonFiniteDirection.into());
    }
    let before = geometry(before_m)?;
    let after = geometry(after_m)?;
    let initial = face.at_area(before.area)?;
    let final_state = face.at_area(after.area)?;
    if face.amount_mol == 0.0 {
        return Ok(triangle_surface_tension_step_tangent(
            before_m,
            after_m,
            initial.tension_n_m,
            initial.tension_n_m,
            direction_m,
        )?);
    }
    let (secant, derivative) = area_secant(face, before.area, after.area, initial, final_state)?;
    let unit = triangle_surface_tension_step(before_m, after_m, 1.0, 1.0)?;
    let mut tangent =
        triangle_surface_tension_step_tangent(before_m, after_m, secant, secant, direction_m)?;
    let da = sub(direction_m[1], direction_m[0]);
    let db = sub(direction_m[2], direction_m[0]);
    let darea =
        0.5 * (dot(cross(da, after.b), after.normal) + dot(cross(after.a, db), after.normal));
    let coefficient_direction = derivative * darea;
    for (force, area_force) in tangent.iter_mut().zip(unit.forces_n) {
        for d in 0..3 {
            force[d] += coefficient_direction * area_force[d];
        }
    }
    if tangent.iter().flatten().any(|v| !v.is_finite()) {
        return Err(SurfactantError::Unrepresentable);
    }
    Ok(tangent)
}

/// Assemble a complete material patch with deterministic scalar accounting.
pub fn assemble_surfactant_step(
    before_m: &[Vector],
    after_m: &[Vector],
    faces: &[[usize; 3]],
    surfactant: SurfactantField<'_>,
) -> Result<CapillaryPatchResponse, SurfactantError> {
    if before_m.len() != after_m.len() {
        return Err(CapillarySurfaceError::LengthMismatch.into());
    }
    surfactant.validate(faces.len())?;
    if before_m
        .iter()
        .chain(after_m)
        .flatten()
        .any(|v| !v.is_finite())
    {
        return Err(CapillarySurfaceError::NonFinitePosition.into());
    }
    let mut forces = vec![[0.0; 3]; before_m.len()];
    let mut totals = [Sum::default(); 7];
    for (face, &indices) in faces.iter().enumerate() {
        for vertex in indices {
            if vertex >= before_m.len() {
                return Err(CapillarySurfaceError::InvalidVertexIndex { face, vertex }.into());
            }
        }
        let before = indices.map(|i| before_m[i]);
        let after = indices.map(|i| after_m[i]);
        let inventory = surfactant.at(face);
        let local = if inventory.amount_mol == 0.0 && before == after {
            triangle_surface_tension(before, inventory.law.parameters().0)?
        } else {
            triangle_surfactant_step(before, after, inventory)?
        };
        for (node, force) in indices.into_iter().zip(local.forces_n) {
            for (total, value) in forces[node].iter_mut().zip(force) {
                *total += value;
                if !total.is_finite() {
                    return Err(SurfactantError::Unrepresentable);
                }
            }
        }
        for (total, value) in totals.iter_mut().zip(local.energy.values()) {
            total.add(value);
        }
    }
    let values = totals.map(Sum::value);
    if values.iter().any(|v| !v.is_finite()) {
        return Err(SurfactantError::Unrepresentable);
    }
    Ok(CapillaryPatchResponse {
        forces_n: forces,
        energy: SurfaceEnergyChange::from_values(values),
    })
}

fn amount(n: f64) -> Result<(), SurfactantError> {
    if !n.is_finite() || n < 0.0 {
        Err(SurfactantError::InvalidAmount)
    } else {
        Ok(())
    }
}

// Compensate the rounded addition before calling the existing deterministic
// logarithm. Its first-order correction is accurate to floating roundoff.
fn log1p(x: f64) -> f64 {
    let y = 1.0 + x;
    if y == 1.0 {
        x
    } else {
        det::ln(y) - ((y - 1.0) - x) / y
    }
}

fn log_ratio(a: f64, b: f64) -> f64 {
    let relative = (a - b) / b;
    if relative.abs() < 0.5 {
        log1p(relative)
    } else {
        let ratio = a / b;
        if ratio.is_finite() && ratio > 0.0 {
            det::ln(ratio)
        } else {
            det::ln(a) - det::ln(b)
        }
    }
}

// Analytic series at removable singularities; at |x|<=1/8 the first omitted
// term is below 1e-24. These are evaluations of the exact EOS secants, not
// finite differences or quadrature replacing the force/tangent law.
fn log1p_ratio(x: f64) -> f64 {
    if x.abs() <= 0.125 {
        let mut value = 0.0;
        for n in (0..=26).rev() {
            value = (-x).mul_add(value, 1.0 / (n + 1) as f64);
        }
        value
    } else {
        log1p(x) / x
    }
}

// H(x)=((1+x)*ln(1+x)-x)/x, with H(0)=0 and H(-1)=-1.
fn h(x: f64) -> f64 {
    if x == -1.0 {
        -1.0
    } else if x.abs() <= 0.125 {
        let mut value = 0.0;
        for n in (1..=26).rev() {
            value = (-x).mul_add(value, 1.0 / (n * (n + 1)) as f64);
        }
        x * value
    } else {
        (1.0 + 1.0 / x) * log1p(x) - 1.0
    }
}

fn h_prime(x: f64) -> f64 {
    if x.abs() <= 0.125 {
        let mut value = 0.0;
        for n in (0..=26).rev() {
            value = (-x).mul_add(value, 1.0 / (n + 2) as f64);
        }
        value
    } else {
        (1.0 - log1p(x) / x) / x
    }
}

// L=ln(A1/A0)/(A1-A0), K=-dL/dA1. Both remain finite at equal areas.
fn logarithmic_area_slope(a0: f64, a1: f64) -> (f64, f64) {
    let change = a1 - a0;
    let x = change / a0;
    if x.abs() <= 0.125 {
        let mut derivative = 0.0;
        for n in (0..=26).rev() {
            derivative = (-x).mul_add(derivative, (n + 1) as f64 / (n + 2) as f64);
        }
        (log1p_ratio(x) / a0, derivative / a0 / a0)
    } else {
        let slope = log_ratio(a1, a0) / change;
        (slope, (slope - 1.0 / a1) / change)
    }
}

fn area_secant(
    face: SurfactantFace,
    a0: f64,
    a1: f64,
    initial: SurfactantState,
    final_state: SurfactantState,
) -> Result<(f64, f64), SurfactantError> {
    if a0 == a1 {
        return Ok((
            initial.tension_n_m,
            0.5 * initial.dilatational_elasticity_n_m / a0,
        ));
    }
    let (clean, rt, reference) = face.law.parameters();
    let amount_energy = rt * face.amount_mol;
    let (slope, derivative) = logarithmic_area_slope(a0, a1);
    let (secant, tangent) = match face.law {
        SurfactantLaw::IdealDilute { .. } => {
            (clean - amount_energy * slope, amount_energy * derivative)
        }
        SurfactantLaw::Langmuir { .. } => {
            let theta0 = initial.concentration_mol_m2 / reference;
            let theta1 = final_state.concentration_mol_m2 / reference;
            if theta0.max(theta1) <= 0.125 {
                // Factor N explicitly so the dilute limit survives even when
                // both packing fractions lie below one ulp of unity.
                let mut s = 1.0;
                let mut t = 1.0;
                let mut power = 1.0;
                let mut sum = Sum::default();
                let mut tangent_sum = Sum::default();
                for n in 0..=26 {
                    let denominator = ((n + 1) * (n + 2)) as f64;
                    sum.add(s / denominator);
                    tangent_sum.add(t / denominator);
                    power *= theta1;
                    s = theta0.mul_add(s, power);
                    t = theta0.mul_add(t, (n + 2) as f64 * power);
                }
                (
                    clean - amount_energy * (slope + theta0 / a1 * sum.value()),
                    amount_energy * (derivative + theta0 / a1 / a1 * tangent_sum.value()),
                )
            } else {
                let change = a1 - a0;
                let x = change / a0;
                let y = change / (a0 * (1.0 - theta0));
                let osmotic = rt * initial.concentration_mol_m2;
                if !x.is_finite() || !y.is_finite() || x <= -1.0 || y <= -1.0 {
                    // Widely separated areas cannot suffer subtractive energy
                    // cancellation; retain their representable endpoint law.
                    let sigma = (final_state.free_energy_j - initial.free_energy_j) / change;
                    (sigma, (final_state.tension_n_m - sigma) / change)
                } else {
                    (
                        clean + (osmotic / theta0) * (log1p(-theta0) + h(y) - h(x)),
                        (osmotic / theta0) * (h_prime(y) / (1.0 - theta0) - h_prime(x)) / a0,
                    )
                }
            }
        }
    };
    if !secant.is_finite() || !tangent.is_finite() || tangent < 0.0 || secant < 0.0 {
        return Err(SurfactantError::Unrepresentable);
    }
    Ok((secant, tangent))
}

#[cfg(test)]
mod tests;
