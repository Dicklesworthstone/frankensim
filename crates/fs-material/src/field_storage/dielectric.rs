//! Passive multipole Debye polarization from an explicitly selected material.
//!
//! With the `exp(i omega t)` convention,
//! `epsilon(omega) = epsilon_infinity + sum(delta_epsilon_k/(1+i omega tau_k))`.
//! A uniform gap therefore has C_infinity in parallel with series R_k/C_k
//! branches, where `C_k = delta_epsilon_k A/d` and `R_k = tau_k/C_k`.
//! Their resistor loss is polarization heat, not DC leakage or electrode ESR.
//! The native circuit owns ALL charges, time integration and dissipation.
//!
//! Sources: COMSOL AC/DC 6.3, "Dispersion Theory", Debye and multipole model;
//! <https://doc.comsol.com/6.3/doc/com.comsol.help.acdc/acdc_ug_theory.05.36.html>.
//! Analog Devices AN-348, Figure 1, gives the series-RC absorption branch.
//! No pole fitting, empirical defaults or measured-accuracy claim is supplied.

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::Time;

use super::{ELECTRIC_PERMITTIVITY_DIMS, FieldStorageError, coefficient, energy, positive};
use crate::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};

/// Source-pinned temperature dependence and complete dielectric thermodynamics.
pub mod thermal;

/// Explicit bound on retained polarization states and per-response work.
pub const MAX_DIELECTRIC_POLES: usize = 16;

/// Original material property names for one polarization process. An increment
/// is absolute permittivity [F/m], not relative epsilon or total static epsilon.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DebyePoleProperties {
    /// Positive absolute polarization increment [F/m].
    pub permittivity_increment: String,
    /// Positive relaxation time [s].
    pub relaxation_time: String,
}

/// One source-derived native series RC branch. These are frozen coefficients at
/// the retained complete material state; changing them needs parameter work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DebyePole {
    /// Original source increment, before geometric reduction [F/m].
    pub permittivity_increment_f_m: f64,
    /// Original source relaxation time [s].
    pub relaxation_time_s: f64,
    /// Polarization storage `delta_epsilon A/d` [F].
    pub capacitance_f: f64,
    /// Polarization friction `tau/C_k` [ohm], distinct from bulk leakage.
    pub resistance_ohm: f64,
}

/// Terminal admittance in siemens, with positive loss for exp(i omega t).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DebyeAdmittance {
    /// Nonnegative real admittance [S].
    pub conductance_s: f64,
    /// Nonnegative imaginary admittance [S].
    pub susceptance_s: f64,
}

/// Immutable uniform-gap geometry and the complete original material bundle.
/// All parameters are explicitly named positive source properties at fixed
/// external conditions. Frequency-, field-, strain- or stress-dependent source
/// parameters cannot be silently frozen into this linear time-domain model.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedDebyeCapacitor {
    material: ResolvedMaterialStatePoint,
    high_frequency_property: String,
    pole_properties: Vec<DebyePoleProperties>,
    area_m2: f64,
    gap_m: f64,
    instantaneous_capacitance_f: f64,
    equilibrium_capacitance_f: f64,
    poles: Vec<DebyePole>,
}

fn invalid(message: &'static str) -> FieldStorageError {
    FieldStorageError::Invalid(message)
}
fn finite(value: f64) -> Result<f64, FieldStorageError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("dielectric response is unrepresentable"))
    }
}
fn admit_names(high: &str, poles: &[DebyePoleProperties]) -> Result<(), FieldStorageError> {
    if high.is_empty() || poles.len() > MAX_DIELECTRIC_POLES {
        return Err(invalid(
            "a named high-frequency source and bounded dielectric pole count are required",
        ));
    }
    let mut names = vec![high];
    for pole in poles {
        for name in [&pole.permittivity_increment, &pole.relaxation_time] {
            if name.is_empty() || names.contains(&name.as_str()) {
                return Err(invalid(
                    "each dielectric parameter needs a distinct original source property",
                ));
            }
            names.push(name);
        }
    }
    Ok(())
}

impl ResolvedDebyeCapacitor {
    /// Adopt an existing complete bundle without dropping other physics claims.
    /// An empty pole list is the lossless high-frequency capacitor. Missing,
    /// zero or negative pole parameters refuse instead of being fitted/clipped.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        high_frequency_property: &str,
        pole_properties: &[DebyePoleProperties],
        area_m2: f64,
        gap_m: f64,
    ) -> Result<Self, FieldStorageError> {
        admit_names(high_frequency_property, pole_properties)?;
        positive(area_m2, "dielectric area must be finite and positive")?;
        positive(gap_m, "dielectric gap must be finite and positive")?;
        let scale = positive(
            area_m2 / gap_m,
            "dielectric geometry ratio is unrepresentable",
        )?;
        let epsilon = coefficient(
            &material,
            high_frequency_property,
            ELECTRIC_PERMITTIVITY_DIMS,
        )?;
        let capacitance = |epsilon| -> Result<f64, FieldStorageError> {
            let c = positive(epsilon * scale, "dielectric capacitance is unrepresentable")?;
            positive(
                c.recip(),
                "inverse dielectric capacitance is unrepresentable",
            )?;
            Ok(c)
        };
        let instantaneous_capacitance_f = capacitance(epsilon)?;
        let mut equilibrium_capacitance_f = instantaneous_capacitance_f;
        let mut poles = Vec::with_capacity(pole_properties.len());
        for properties in pole_properties {
            let delta = coefficient(
                &material,
                &properties.permittivity_increment,
                ELECTRIC_PERMITTIVITY_DIMS,
            )?;
            let relaxation_time_s =
                coefficient(&material, &properties.relaxation_time, Time::DIMS)?;
            let capacitance_f = capacitance(delta)?;
            let resistance_ohm = positive(
                relaxation_time_s / capacitance_f,
                "dielectric relaxation resistance is unrepresentable",
            )?;
            positive(
                resistance_ohm.recip(),
                "dielectric relaxation conductance is unrepresentable",
            )?;
            equilibrium_capacitance_f = finite(equilibrium_capacitance_f + capacitance_f)?;
            poles.push(DebyePole {
                permittivity_increment_f_m: delta,
                relaxation_time_s,
                capacitance_f,
                resistance_ohm,
            });
        }
        Ok(Self {
            material,
            high_frequency_property: high_frequency_property.to_owned(),
            pole_properties: pole_properties.to_vec(),
            area_m2,
            gap_m,
            instantaneous_capacitance_f,
            equilibrium_capacitance_f,
            poles,
        })
    }

    /// Resolve exactly the declared parameter roles through the common material
    /// query owner. Source identities, validity and physical coordinate meanings
    /// remain attached; the caller explicitly declares the Debye model semantics.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        high_frequency_property: &str,
        pole_properties: &[DebyePoleProperties],
        area_m2: f64,
        gap_m: f64,
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
            area_m2,
            gap_m,
        )
    }

    /// Re-evaluate the entire original material bundle at a new query, retaining
    /// every selected claim, coordinate convention, pole role and fixed geometry.
    /// This supplies an endpoint law; it does not change any electric history or
    /// account for the work/entropy of changing the coefficients.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, FieldStorageError> {
        Self::try_from_material(
            self.material.at_query_point(card, point)?,
            &self.high_frequency_property,
            &self.pole_properties,
            self.area_m2,
            self.gap_m,
        )
    }

    /// Temperature-parameter change in Helmholtz energy at fixed instantaneous
    /// and polarization charges. Geometry, the complete source bundle and all
    /// nonthermal coordinates must remain identical. This signed quantity is
    /// reversible parameter work, not irreversible polarization heat or total
    /// internal-energy change; thermal owners also need the entropy correction.
    ///
    /// Evaluate each contribution as the larger endpoint energy times the
    /// relative capacitance change, avoiding subtraction of large energy values
    /// or nearly equal reciprocals. Opposite pole contributions are compensated.
    pub fn fixed_charge_free_energy_change_j(
        &self,
        refreshed: &Self,
        instantaneous_charge_c: f64,
        polarization_charges_c: &[f64],
    ) -> Result<f64, FieldStorageError> {
        super::admit_source_temperature_refresh(&self.material, &refreshed.material)?;
        if self.high_frequency_property != refreshed.high_frequency_property
            || self.pole_properties != refreshed.pole_properties
            || self.area_m2.to_bits() != refreshed.area_m2.to_bits()
            || self.gap_m.to_bits() != refreshed.gap_m.to_bits()
        {
            return Err(invalid(
                "dielectric parameter work requires the same pole roles and fixed geometry",
            ));
        }
        self.stored_energy_j(instantaneous_charge_c, polarization_charges_c)?;
        refreshed.stored_energy_j(instantaneous_charge_c, polarization_charges_c)?;
        let mut sum = 0.0_f64;
        let mut correction = 0.0_f64;
        let mut add = |charge, before: f64, after: f64| -> Result<(), FieldStorageError> {
            let term = if before == after || charge == 0.0 {
                0.0
            } else {
                let larger_energy = energy(charge, before.min(after))?;
                let change = finite(larger_energy * ((before - after) / before.max(after)))?;
                if change == 0.0 {
                    return Err(invalid("nonzero dielectric parameter work underflows"));
                }
                change
            };
            let next = finite(sum + term)?;
            correction = finite(
                correction
                    + if sum.abs() >= term.abs() {
                        (sum - next) + term
                    } else {
                        (term - next) + sum
                    },
            )?;
            sum = next;
            Ok(())
        };
        add(
            instantaneous_charge_c,
            self.instantaneous_capacitance_f,
            refreshed.instantaneous_capacitance_f,
        )?;
        for ((before, after), &charge) in self
            .poles
            .iter()
            .zip(&refreshed.poles)
            .zip(polarization_charges_c)
        {
            add(charge, before.capacitance_f, after.capacitance_f)?;
        }
        finite(sum + correction)
    }

    /// Complete original source bundle and physical query.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Explicit original property selected for epsilon_infinity.
    pub fn high_frequency_property(&self) -> &str {
        &self.high_frequency_property
    }
    /// Original parameter-property names, in retained polarization order.
    pub fn pole_properties(&self) -> &[DebyePoleProperties] {
        &self.pole_properties
    }
    /// Fixed facing area [m2].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Fixed dielectric thickness [m].
    pub const fn gap_m(&self) -> f64 {
        self.gap_m
    }
    /// High-frequency capacitance with all polarization processes frozen [F].
    pub const fn instantaneous_capacitance_f(&self) -> f64 {
        self.instantaneous_capacitance_f
    }
    /// Static capacitance after all polarization processes equilibrate [F].
    pub const fn equilibrium_capacitance_f(&self) -> f64 {
        self.equilibrium_capacitance_f
    }
    /// Exact source-derived RC coefficients, in retained polarization order.
    pub fn poles(&self) -> &[DebyePole] {
        &self.poles
    }

    /// Allocation-free `Y = i omega C_inf + sum(i omega C_k/(1+i omega tau_k))`.
    /// Angular frequency is in rad/s. Stable reciprocal evaluation preserves the
    /// high-frequency conductance limit even when omega*tau overflows. DC loss
    /// is exactly zero: polarization relaxation does not create leakage.
    pub fn admittance(&self, omega_rad_s: f64) -> Result<DebyeAdmittance, FieldStorageError> {
        if !omega_rad_s.is_finite() || omega_rad_s < 0.0 {
            return Err(invalid(
                "dielectric admittance needs finite nonnegative angular frequency",
            ));
        }
        let mut conductance_s = 0.0;
        let mut susceptance_s = finite(omega_rad_s * self.instantaneous_capacitance_f)?;
        for pole in &self.poles {
            let x = omega_rad_s * pole.relaxation_time_s;
            let (g, b) = if x <= 1.0 {
                let b = omega_rad_s * pole.capacitance_f / (1.0 + x * x);
                (b * x, b)
            } else {
                let inverse = (1.0 / omega_rad_s) / pole.relaxation_time_s;
                let g = pole.resistance_ohm.recip() / (1.0 + inverse * inverse);
                (g, g * inverse)
            };
            conductance_s = finite(conductance_s + g)?;
            susceptance_s = finite(susceptance_s + b)?;
        }
        Ok(DebyeAdmittance {
            conductance_s,
            susceptance_s,
        })
    }

    /// Physical electric storage for the instantaneous charge and each retained
    /// polarization charge, in pole order. Terminal charge is their SUM; using
    /// q_total squared over a frequency-dependent C would erase the memory.
    pub fn stored_energy_j(
        &self,
        instantaneous_charge_c: f64,
        polarization_charges_c: &[f64],
    ) -> Result<f64, FieldStorageError> {
        if polarization_charges_c.len() != self.poles.len() {
            return Err(invalid(
                "one original charge per dielectric polarization pole required",
            ));
        }
        let mut result = energy(instantaneous_charge_c, self.instantaneous_capacitance_f)?;
        for (pole, &charge) in self.poles.iter().zip(polarization_charges_c) {
            result = finite(result + energy(charge, pole.capacitance_f)?)?;
        }
        Ok(result)
    }
}
