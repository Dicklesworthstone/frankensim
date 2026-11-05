//! Local multipole Debye memory for spatial electroquasistatic fields.
//!
//! Each fixed isotropic pole obeys `tau Pdot + P = delta_epsilon E`, with
//! `D = epsilon_infinity E + sum P`. The free-energy density is
//! `epsilon_infinity |E|²/2 + sum |P|²/(2 delta_epsilon)`, and the irreversible
//! heat rate is `sum tau |Pdot|²/delta_epsilon`. Polarization is not an Ohmic
//! current and must not also be put into the conductivity coefficient.
//!
//! Midpoint integration eliminates each local pole before the field solve.
//! With `b=dt/(2 tau+dt)`, the endpoint relation is
//! `D1=(epsilon_infinity+sum b delta_epsilon) E1 + history`, where
//! `history=sum [(1-2b)P0+b delta_epsilon E0]`. Recovery uses the same b and
//! exposes storage, work and nonnegative polarization heat independently.
//! This is second-order, passive midpoint integration, not an exact exponential
//! step or a sign-preserving scheme at arbitrary dt; callers must resolve the
//! relevant relaxation times and independently refine space and time.
//!
//! Constitutive reference: COMSOL 6.3, "Dispersion Theory", multipole Debye and
//! local polarization variables:
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.36.html>.
//! The storage/work identity above follows by multiplying each pole equation
//! by `Pdot/delta_epsilon`; no fit or measured material accuracy is implied.

use super::{
    DebyePoleProperties, ELECTRIC_PERMITTIVITY_DIMS, FieldStorageError, MaterialCard,
    MaterialPropertySelection, QueryPoint, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, Time, admit_names, coefficient, energy, finite, invalid,
    resolve_material_state_point,
};

/// One original bulk polarization process, before any geometric reduction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkDebyePole {
    /// Positive absolute increment in permittivity [F/m], not relative epsilon.
    pub permittivity_increment_f_m: f64,
    /// Positive relaxation time [s].
    pub relaxation_time_s: f64,
}

/// Fixed-condition source-backed bulk dielectric; all coefficients and the
/// complete original material bundle stay immutable throughout a trajectory.
/// The named roles share the lumped Debye source/quantity/domain admission.
/// There is no synthetic unit capacitor, inferred pole, or parameter refresh.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedDebyeMaterial {
    material: ResolvedMaterialStatePoint,
    high_frequency_property: String,
    pole_properties: Vec<DebyePoleProperties>,
    high_frequency_permittivity_f_m: f64,
    equilibrium_permittivity_f_m: f64,
    poles: Vec<BulkDebyePole>,
}

/// Local recovery for one accepted field interval. Energy quantities are
/// densities [J/m3]; the spatial owner integrates them over original volumes.
#[derive(Clone, Debug, PartialEq)]
pub struct DebyeMidpointStep<const N: usize> {
    /// Endpoint polarization, in the original material pole order [C/m2].
    pub polarization_after_c_m2: Vec<[f64; N]>,
    /// Total original displacement, including every polarization pole [C/m2].
    pub displacement_before_c_m2: [f64; N],
    /// Total endpoint displacement, including every polarization pole [C/m2].
    pub displacement_after_c_m2: [f64; N],
    /// Original physical storage density [J/m3].
    pub stored_before_j_m3: f64,
    /// Endpoint physical storage density [J/m3].
    pub stored_after_j_m3: f64,
    /// Stable change in storage computed from actual endpoint increments.
    pub stored_change_j_m3: f64,
    /// Nonnegative irreversible polarization heat, excluding Ohmic heating.
    pub heat_j_m3: f64,
    /// Midpoint electric work `E_mid . (D1-D0)`, using component increments.
    pub supplied_work_j_m3: f64,
    /// `stored_change + heat - supplied_work`; callers own physical tolerances.
    pub energy_residual_j_m3: f64,
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), FieldStorageError> {
        let next = finite(self.value + finite(term)?)?;
        self.correction = finite(
            self.correction
                + if self.value.abs() >= term.abs() {
                    (self.value - next) + term
                } else {
                    (term - next) + self.value
                },
        )?;
        self.value = next;
        Ok(())
    }
    fn total(&self) -> Result<f64, FieldStorageError> {
        finite(self.value + self.correction)
    }
}

fn interval(dt_s: f64) -> Result<(), FieldStorageError> {
    if dt_s.is_finite() && dt_s > 0.0 {
        Ok(())
    } else {
        Err(invalid(
            "bulk dielectric interval must be finite and positive",
        ))
    }
}

// Both b and its complement are evaluated directly. Computing 1-b would erase
// the small positive loss of a pole much faster than the chosen interval.
fn weights(dt_s: f64, tau_s: f64) -> (f64, f64) {
    if dt_s <= tau_s {
        let ratio = dt_s / tau_s;
        let denominator = 2.0 + ratio;
        (ratio / denominator, 2.0 / denominator)
    } else {
        let ratio = tau_s / dt_s;
        let denominator = 1.0 + 2.0 * ratio;
        (denominator.recip(), 2.0 * ratio / denominator)
    }
}

fn field_energy(e: f64, epsilon: f64) -> Result<f64, FieldStorageError> {
    let displacement = finite(epsilon * e)?;
    let value = finite(0.5 * e * displacement)?;
    if e != 0.0 && (displacement == 0.0 || value == 0.0) {
        return Err(invalid("nonzero bulk dielectric storage underflows"));
    }
    Ok(value)
}

impl ResolvedDebyeMaterial {
    /// Adopt the complete common material bundle with explicit epsilon_infinity
    /// and pole roles. An empty pole list is exactly the lossless linear law.
    /// Positive absolute SI coefficients and fixed external conditions are
    /// mandatory; frequency/field/strain/stress-dependent parameters refuse.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        high_frequency_property: &str,
        pole_properties: &[DebyePoleProperties],
    ) -> Result<Self, FieldStorageError> {
        admit_names(high_frequency_property, pole_properties)?;
        let high_frequency_permittivity_f_m = coefficient(
            &material,
            high_frequency_property,
            ELECTRIC_PERMITTIVITY_DIMS,
        )?;
        let mut equilibrium = Sum::default();
        equilibrium.add(high_frequency_permittivity_f_m)?;
        let mut poles = Vec::with_capacity(pole_properties.len());
        for names in pole_properties {
            let permittivity_increment_f_m = coefficient(
                &material,
                &names.permittivity_increment,
                ELECTRIC_PERMITTIVITY_DIMS,
            )?;
            let relaxation_time_s = coefficient(&material, &names.relaxation_time, Time::DIMS)?;
            equilibrium.add(permittivity_increment_f_m)?;
            poles.push(BulkDebyePole {
                permittivity_increment_f_m,
                relaxation_time_s,
            });
        }
        Ok(Self {
            material,
            high_frequency_property: high_frequency_property.to_owned(),
            pole_properties: pole_properties.to_vec(),
            high_frequency_permittivity_f_m,
            equilibrium_permittivity_f_m: equilibrium.total()?,
            poles,
        })
    }

    /// Resolve each original parameter through the shared source owner without
    /// fitting poles or adding material properties at the field boundary.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        high_frequency_property: &str,
        pole_properties: &[DebyePoleProperties],
    ) -> Result<Self, FieldStorageError> {
        admit_names(high_frequency_property, pole_properties)?;
        let requirement = |name: &str, dims| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
        };
        let mut requirements = vec![requirement(
            high_frequency_property,
            ELECTRIC_PERMITTIVITY_DIMS,
        )?];
        for pole in pole_properties {
            requirements.push(requirement(
                &pole.permittivity_increment,
                ELECTRIC_PERMITTIVITY_DIMS,
            )?);
            requirements.push(requirement(&pole.relaxation_time, Time::DIMS)?);
        }
        Self::try_from_material(
            resolve_material_state_point(card, point, &requirements, selection)?,
            high_frequency_property,
            pole_properties,
        )
    }

    /// Complete original source bundle at its fixed physical query.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Original named source for epsilon_infinity.
    pub fn high_frequency_property(&self) -> &str {
        &self.high_frequency_property
    }
    /// Original named parameter roles, in polarization-state order.
    pub fn pole_properties(&self) -> &[DebyePoleProperties] {
        &self.pole_properties
    }
    /// Absolute instantaneous bulk permittivity [F/m].
    pub const fn high_frequency_permittivity_f_m(&self) -> f64 {
        self.high_frequency_permittivity_f_m
    }
    /// Absolute permittivity after all poles equilibrate [F/m].
    pub const fn equilibrium_permittivity_f_m(&self) -> f64 {
        self.equilibrium_permittivity_f_m
    }
    /// Original bulk coefficients, with no capacitor geometry attached.
    pub fn poles(&self) -> &[BulkDebyePole] {
        &self.poles
    }

    fn admit_state<const N: usize>(
        &self,
        field: &[f64; N],
        polarization: &[[f64; N]],
    ) -> Result<(), FieldStorageError> {
        if !(1..=3).contains(&N)
            || polarization.len() != self.poles.len()
            || field
                .iter()
                .chain(polarization.iter().flatten())
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "bulk dielectric requires one finite 1D/2D/3D polarization vector per original pole",
            ));
        }
        Ok(())
    }

    /// Total displacement at an explicit local field/history, without allocation.
    pub fn displacement_c_m2<const N: usize>(
        &self,
        electric_v_m: [f64; N],
        polarization_c_m2: &[[f64; N]],
    ) -> Result<[f64; N], FieldStorageError> {
        self.admit_state(&electric_v_m, polarization_c_m2)?;
        let mut result = [0.0; N];
        for i in 0..N {
            let mut sum = Sum::default();
            sum.add(self.high_frequency_permittivity_f_m * electric_v_m[i])?;
            for p in polarization_c_m2 {
                sum.add(p[i])?;
            }
            result[i] = sum.total()?;
        }
        Ok(result)
    }

    /// Physical storage `epsilon_inf |E|²/2 + sum |P_k|²/(2 delta_epsilon_k)`.
    /// It depends on each pole separately, not on total D and a fitted epsilon.
    pub fn stored_energy_density_j_m3<const N: usize>(
        &self,
        electric_v_m: [f64; N],
        polarization_c_m2: &[[f64; N]],
    ) -> Result<f64, FieldStorageError> {
        self.admit_state(&electric_v_m, polarization_c_m2)?;
        let mut sum = Sum::default();
        for e in electric_v_m {
            sum.add(field_energy(e, self.high_frequency_permittivity_f_m)?)?;
        }
        for (pole, p) in self.poles.iter().zip(polarization_c_m2) {
            for value in p {
                sum.add(energy(*value, pole.permittivity_increment_f_m)?)?;
            }
        }
        sum.total()
    }

    /// Positive coefficient for the endpoint field solve after local midpoint
    /// elimination. O(poles), no allocation or repeated material query.
    pub fn midpoint_permittivity_f_m(&self, dt_s: f64) -> Result<f64, FieldStorageError> {
        interval(dt_s)?;
        let mut sum = Sum::default();
        sum.add(self.high_frequency_permittivity_f_m)?;
        for pole in &self.poles {
            let (b, _) = weights(dt_s, pole.relaxation_time_s);
            sum.add(b * pole.permittivity_increment_f_m)?;
        }
        sum.total()
    }

    /// Allocation-free history in `D1 = epsilon_effective E1 + history`.
    /// Use this SAME material and interval for endpoint recovery.
    pub fn midpoint_history_c_m2<const N: usize>(
        &self,
        dt_s: f64,
        electric_before_v_m: [f64; N],
        polarization_before_c_m2: &[[f64; N]],
    ) -> Result<[f64; N], FieldStorageError> {
        interval(dt_s)?;
        self.admit_state(&electric_before_v_m, polarization_before_c_m2)?;
        let mut result = [0.0; N];
        for i in 0..N {
            let mut sum = Sum::default();
            for (pole, p) in self.poles.iter().zip(polarization_before_c_m2) {
                let (b, complement) = weights(dt_s, pole.relaxation_time_s);
                sum.add((complement - b).mul_add(
                    p[i],
                    (b * pole.permittivity_increment_f_m) * electric_before_v_m[i],
                ))?;
            }
            result[i] = sum.total()?;
        }
        Ok(result)
    }

    /// Recover all local polarization endpoints after the coupled field solve.
    /// Only the at-most-16-pole output vector allocates. The input history is
    /// untouched on failure; geometry integration and tolerance admission belong
    /// to the caller. A small energy residual includes rounding of retained P.
    pub fn step_midpoint<const N: usize>(
        &self,
        dt_s: f64,
        electric_before_v_m: [f64; N],
        electric_after_v_m: [f64; N],
        polarization_before_c_m2: &[[f64; N]],
    ) -> Result<DebyeMidpointStep<N>, FieldStorageError> {
        interval(dt_s)?;
        self.admit_state(&electric_before_v_m, polarization_before_c_m2)?;
        self.admit_state(&electric_after_v_m, polarization_before_c_m2)?;
        let midpoint: [f64; N] =
            std::array::from_fn(|i| f64::midpoint(electric_before_v_m[i], electric_after_v_m[i]));
        let mut polarization_after_c_m2 = Vec::with_capacity(self.poles.len());
        let mut change = Sum::default();
        let mut work = Sum::default();
        let mut heat = Sum::default();
        for i in 0..N {
            let delta = finite(electric_after_v_m[i] - electric_before_v_m[i])?;
            let term = finite((self.high_frequency_permittivity_f_m * delta) * midpoint[i])?;
            change.add(term)?;
            work.add(term)?;
        }
        for (pole, before) in self.poles.iter().zip(polarization_before_c_m2) {
            let (b, complement) = weights(dt_s, pole.relaxation_time_s);
            let epsilon = pole.permittivity_increment_f_m;
            let mut after = [0.0; N];
            for i in 0..N {
                let gap = finite(epsilon.mul_add(midpoint[i], -before[i]))?;
                after[i] = finite((2.0 * b).mul_add(gap, before[i]))?;
                let increment = finite(after[i] - before[i])?;
                change.add((increment / epsilon) * f64::midpoint(before[i], after[i]))?;
                work.add(midpoint[i] * increment)?;
                // 4 b(1-b) * gap²/(2 epsilon) = tau |Delta P|²/(epsilon dt)
                // for the unrounded midpoint update. Direct complement retains
                // positive heat when b rounds to one in a very stiff interval.
                heat.add((4.0 * b * complement) * energy(gap, epsilon)?)?;
            }
            polarization_after_c_m2.push(after);
        }
        let stored_change_j_m3 = change.total()?;
        let heat_j_m3 = heat.total()?;
        let supplied_work_j_m3 = work.total()?;
        let mut residual = Sum::default();
        residual.add(stored_change_j_m3)?;
        residual.add(heat_j_m3)?;
        residual.add(-supplied_work_j_m3)?;
        Ok(DebyeMidpointStep {
            displacement_before_c_m2: self
                .displacement_c_m2(electric_before_v_m, polarization_before_c_m2)?,
            displacement_after_c_m2: self
                .displacement_c_m2(electric_after_v_m, &polarization_after_c_m2)?,
            stored_before_j_m3: self
                .stored_energy_density_j_m3(electric_before_v_m, polarization_before_c_m2)?,
            stored_after_j_m3: self
                .stored_energy_density_j_m3(electric_after_v_m, &polarization_after_c_m2)?,
            polarization_after_c_m2,
            stored_change_j_m3,
            heat_j_m3,
            supplied_work_j_m3,
            energy_residual_j_m3: residual.total()?,
        })
    }
}
