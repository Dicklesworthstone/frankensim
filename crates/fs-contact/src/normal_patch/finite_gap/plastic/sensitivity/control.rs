//! Fixed-history tangents under physical force and finite-spring boundaries.

use fs_exec::Cx;

use super::super::{
    FiniteGapPlasticDrive, FiniteGapPlasticDriveAdvance, FiniteGapPlasticLoadAdvance,
    load::validate_controls,
};
use super::{
    FiniteGapPlasticDirection, FiniteGapPlasticLinearization, FiniteGapPlasticModel,
    FiniteGapPlasticSensitivity, FiniteGapPlasticSensitivityControls,
    FiniteGapPlasticSensitivityError, FiniteGapPlasticSensitivityMargins, FiniteGapPlasticState,
    checkpoint,
};

/// Accuracy and local uniqueness required for a prescribed-force tangent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGapPlasticLoadSensitivityControls {
    /// Active-set separation and physical free-equation accuracy.
    pub contact: FiniteGapPlasticSensitivityControls,
    /// Strict positive lower bound on contact force per unit approach [N/m].
    /// A zero-stiffness plastic plateau cannot determine indentation at fixed force.
    pub minimum_contact_stiffness_n_per_m: f64,
    /// Maximum differentiated force-balance defect [N per dimensionless parameter].
    pub maximum_force_residual_n: f64,
}

/// Accuracy required for a finite-loading-spring tangent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGapPlasticDriveSensitivityControls {
    /// Active-set separation and physical free-equation accuracy.
    pub contact: FiniteGapPlasticSensitivityControls,
    /// Maximum differentiated force-balance defect [N per dimensionless parameter].
    pub maximum_force_residual_n: f64,
}

/// Signed force and material/geometry derivatives with accepted history fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGapPlasticLoadDirection<'a> {
    /// Prescribed normal-force derivative [N]. Zero maintains the current load.
    pub target_force_n: f64,
    /// Homogeneous reduced-modulus derivative [Pa].
    pub reduced_modulus_pa: f64,
    /// Pressure-hardness derivative at each original cell [Pa].
    pub hardness_pa: &'a [f64],
    /// Relative recession derivative at each original cell [m].
    pub recession_m: &'a [f64],
}

/// Signed loading-spring and material/geometry derivatives with fixed history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGapPlasticDriveDirection<'a> {
    /// Loading-spring endpoint derivative [m].
    pub position_m: f64,
    /// Loading-spring stiffness derivative [N/m].
    pub stiffness_n_per_m: f64,
    /// Homogeneous reduced-modulus derivative [Pa].
    pub reduced_modulus_pa: f64,
    /// Pressure-hardness derivative at each original cell [Pa].
    pub hardness_pa: &'a [f64],
    /// Relative recession derivative at each original cell [m].
    pub recession_m: &'a [f64],
}

/// Prepared prescribed-force tangent with one cached unit-approach response.
#[derive(Debug)]
pub struct FiniteGapPlasticLoadLinearization<'a> {
    contact: FiniteGapPlasticLinearization<'a>,
    unit_pressure_pa: Vec<f64>,
    contact_stiffness_n_per_m: f64,
    maximum_force_residual_n: f64,
}

/// Prepared finite-spring tangent using the original augmented free operator.
#[derive(Debug)]
pub struct FiniteGapPlasticDriveLinearization<'a> {
    contact: FiniteGapPlasticLinearization<'a>,
    boundary: FiniteGapPlasticDrive,
    compression_m: f64,
    rank_one: f64,
    maximum_force_residual_n: f64,
}

/// Local derivative of one prescribed-force contact advance.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticLoadSensitivity {
    /// Pressure, opening, indentation, force, storage and plastic-work derivatives.
    pub contact: FiniteGapPlasticSensitivity,
    /// Rigid-approach derivative required by the prescribed load [m].
    pub approach_m: f64,
    /// Absolute recomputed differentiated force-balance defect [N].
    pub force_residual_n: f64,
}

/// Local derivative of one finite-loading-spring contact advance.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGapPlasticDriveSensitivity {
    /// Pressure, opening, indentation, force, storage and plastic-work derivatives.
    pub contact: FiniteGapPlasticSensitivity,
    /// Rigid-approach derivative required by the loading spring [m].
    pub approach_m: f64,
    /// Signed loading-spring storage derivative [J], not heat or actuator work.
    pub spring_energy_j: f64,
    /// Absolute recomputed differentiated contact-minus-spring force defect [N].
    pub force_residual_n: f64,
}

impl FiniteGapPlasticModel {
    /// Prepare a tangent from an original ordinary prescribed-force advance.
    ///
    /// The original predecessor, pressure cap and plastic increment are admitted
    /// by the same strict active-set checks as [`Self::linearize_advance`]. One
    /// unit-approach solve supplies the contact stiffness and is reused by every
    /// direction. Fully plastic and separated zero-stiffness states are refused:
    /// force alone does not select a unique approach there. Refresh/recession
    /// continuation receipts are not interchangeable with an ordinary advance.
    pub fn linearize_load_advance<'a>(
        &'a self,
        before: &'a FiniteGapPlasticState,
        advance: &'a FiniteGapPlasticLoadAdvance,
        controls: FiniteGapPlasticLoadSensitivityControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLoadLinearization<'a>, FiniteGapPlasticSensitivityError> {
        admit_force_budget(controls.maximum_force_residual_n)?;
        if !controls.minimum_contact_stiffness_n_per_m.is_finite()
            || controls.minimum_contact_stiffness_n_per_m <= 0.0
        {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "finite positive minimum contact stiffness",
            ));
        }
        let contact = self.linearize_advance(before, &advance.contact, controls.contact, cx)?;
        let receipt = &advance.load;
        if receipt.model_identity != self.identity
            || receipt.predecessor_state_identity != before.identity
            || receipt.candidate_state_identity != advance.contact.candidate.identity
            || receipt.response_identity != advance.contact.receipt.identity
        {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "original load receipt and contact identities",
            ));
        }
        let tolerance =
            validate_controls(receipt.target_force_n, receipt.controls).map_err(|_| {
                FiniteGapPlasticSensitivityError::Invalid(
                    "original finite load controls and target",
                )
            })?;
        let force = pressure_force(self, &contact.candidate.pressure_pa, cx)?;
        if (force - receipt.target_force_n).abs() > tolerance {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "candidate satisfies the original prescribed force",
            ));
        }
        let zero = vec![0.0; contact.opening_m.len()];
        let unit = contact.apply(
            FiniteGapPlasticDirection {
                approach_m: 1.0,
                reduced_modulus_pa: 0.0,
                hardness_pa: &zero,
                recession_m: &zero,
            },
            cx,
        )?;
        let stiffness = unit.normal_force_n;
        if stiffness <= controls.minimum_contact_stiffness_n_per_m {
            return Err(FiniteGapPlasticSensitivityError::UnresolvedLoadStiffness {
                stiffness_n_per_m: stiffness,
                minimum_n_per_m: controls.minimum_contact_stiffness_n_per_m,
            });
        }
        Ok(FiniteGapPlasticLoadLinearization {
            contact,
            unit_pressure_pa: unit.pressure_pa,
            contact_stiffness_n_per_m: stiffness,
            maximum_force_residual_n: controls.maximum_force_residual_n,
        })
    }

    /// Prepare a tangent from an original ordinary finite-spring advance.
    ///
    /// Eliminating approach adds `(cell_area / K) 11^T` to the original free
    /// compliance. This remains well posed on a fully plastic plateau without
    /// inventing contact stiffness. Endpoint and stiffness derivatives are both
    /// supported; they do not describe an unresolved actuator time history.
    pub fn linearize_drive_advance<'a>(
        &'a self,
        before: &'a FiniteGapPlasticState,
        advance: &'a FiniteGapPlasticDriveAdvance,
        controls: FiniteGapPlasticDriveSensitivityControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticDriveLinearization<'a>, FiniteGapPlasticSensitivityError> {
        admit_force_budget(controls.maximum_force_residual_n)?;
        let contact = self.linearize_advance(before, &advance.contact, controls.contact, cx)?;
        let receipt = &advance.drive;
        if receipt.model_identity != self.identity
            || receipt.predecessor_state_identity != before.identity
            || receipt.candidate_state_identity != advance.contact.candidate.identity
            || receipt.response_identity != advance.contact.receipt.identity
        {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "original drive receipt and contact identities",
            ));
        }
        validate_controls(0.0, receipt.controls).map_err(|_| {
            FiniteGapPlasticSensitivityError::Invalid("original finite drive controls")
        })?;
        let boundary = receipt.boundary;
        let compression_m = boundary.position_m - contact.candidate.approach_m;
        let force = pressure_force(self, &contact.candidate.pressure_pa, cx)?;
        let spring_force = boundary.stiffness_n_per_m * compression_m;
        let tolerance = receipt.controls.absolute_force_tolerance_n
            + receipt.controls.relative_force_tolerance * force.abs().max(spring_force.abs());
        let area = self.config.grid.cell_width_m * self.config.grid.cell_depth_m;
        let rank_one = area / boundary.stiffness_n_per_m;
        if !boundary.position_m.is_finite()
            || !boundary.stiffness_n_per_m.is_finite()
            || boundary.stiffness_n_per_m <= 0.0
            || !compression_m.is_finite()
            || compression_m < 0.0
            || !spring_force.is_finite()
            || !tolerance.is_finite()
            || !rank_one.is_finite()
            || rank_one <= 0.0
            || (compression_m != 0.0 && spring_force == 0.0)
            || (force == 0.0 && spring_force != 0.0)
            || (force - spring_force).abs() > tolerance
        {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "candidate satisfies the original finite positive loading spring",
            ));
        }
        Ok(FiniteGapPlasticDriveLinearization {
            contact,
            boundary,
            compression_m,
            rank_one,
            maximum_force_residual_n: controls.maximum_force_residual_n,
        })
    }
}

impl FiniteGapPlasticLoadLinearization<'_> {
    /// Active-set separation of the original controlled equilibrium.
    pub const fn margins(&self) -> FiniteGapPlasticSensitivityMargins {
        self.contact.margins()
    }

    /// Computed contact stiffness with original history and active set fixed [N/m].
    pub const fn contact_stiffness_n_per_m(&self) -> f64 {
        self.contact_stiffness_n_per_m
    }

    /// Apply a signed direction with one original free solve and a scalar balance.
    pub fn apply(
        &self,
        direction: FiniteGapPlasticLoadDirection<'_>,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticLoadSensitivity, FiniteGapPlasticSensitivityError> {
        if !direction.target_force_n.is_finite() {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "finite signed force direction",
            ));
        }
        let fixed = FiniteGapPlasticDirection {
            approach_m: 0.0,
            reduced_modulus_pa: direction.reduced_modulus_pa,
            hardness_pa: direction.hardness_pa,
            recession_m: direction.recession_m,
        };
        let mut pressure = self.contact.solve_direction(fixed, 0.0, cx)?;
        let fixed_force = pressure_force(self.contact.model, &pressure, cx)?;
        let approach_m = (direction.target_force_n - fixed_force) / self.contact_stiffness_n_per_m;
        if !approach_m.is_finite() {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "representable force-controlled approach direction",
            ));
        }
        for (i, value) in pressure.iter_mut().enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            *value = approach_m.mul_add(self.unit_pressure_pa[i], *value);
        }
        // Recheck the combined physical direction: a small unit-response error
        // can otherwise be amplified by a large resolved approach derivative.
        let contact = self.contact.finish_direction(
            FiniteGapPlasticDirection {
                approach_m,
                ..fixed
            },
            pressure,
            cx,
        )?;
        let force_residual_n = check_force_residual(
            contact.normal_force_n - direction.target_force_n,
            self.maximum_force_residual_n,
        )?;
        Ok(FiniteGapPlasticLoadSensitivity {
            contact,
            approach_m,
            force_residual_n,
        })
    }
}

impl FiniteGapPlasticDriveLinearization<'_> {
    /// Active-set separation of the original controlled equilibrium.
    pub const fn margins(&self) -> FiniteGapPlasticSensitivityMargins {
        self.contact.margins()
    }

    /// Apply a signed direction with one spring-augmented native free solve.
    pub fn apply(
        &self,
        direction: FiniteGapPlasticDriveDirection<'_>,
        cx: &Cx<'_>,
    ) -> Result<FiniteGapPlasticDriveSensitivity, FiniteGapPlasticSensitivityError> {
        if !direction.position_m.is_finite() || !direction.stiffness_n_per_m.is_finite() {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "finite signed loading-spring direction",
            ));
        }
        let stiffness = self.boundary.stiffness_n_per_m;
        let effective = FiniteGapPlasticDirection {
            approach_m: self.compression_m.mul_add(
                direction.stiffness_n_per_m / stiffness,
                direction.position_m,
            ),
            reduced_modulus_pa: direction.reduced_modulus_pa,
            hardness_pa: direction.hardness_pa,
            recession_m: direction.recession_m,
        };
        let pressure = self.contact.solve_direction(effective, self.rank_one, cx)?;
        let force_direction = pressure_force(self.contact.model, &pressure, cx)?;
        let approach_m = effective.approach_m - force_direction / stiffness;
        if !approach_m.is_finite() {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "representable finite-spring approach direction",
            ));
        }
        let contact = self.contact.finish_direction(
            FiniteGapPlasticDirection {
                approach_m,
                ..effective
            },
            pressure,
            cx,
        )?;
        let spring_force_direction = stiffness.mul_add(
            direction.position_m - approach_m,
            self.compression_m * direction.stiffness_n_per_m,
        );
        let force_residual_n = check_force_residual(
            contact.normal_force_n - spring_force_direction,
            self.maximum_force_residual_n,
        )?;
        let spring_energy_j = self.compression_m
            * (contact.normal_force_n - (0.5 * self.compression_m) * direction.stiffness_n_per_m);
        if !spring_energy_j.is_finite() {
            return Err(FiniteGapPlasticSensitivityError::Invalid(
                "representable loading-spring storage direction",
            ));
        }
        Ok(FiniteGapPlasticDriveSensitivity {
            contact,
            approach_m,
            spring_energy_j,
            force_residual_n,
        })
    }
}

fn admit_force_budget(maximum_n: f64) -> Result<(), FiniteGapPlasticSensitivityError> {
    if !maximum_n.is_finite() || maximum_n <= 0.0 {
        return Err(FiniteGapPlasticSensitivityError::Invalid(
            "finite positive directional force residual budget",
        ));
    }
    Ok(())
}

fn check_force_residual(
    residual_n: f64,
    maximum_n: f64,
) -> Result<f64, FiniteGapPlasticSensitivityError> {
    if !residual_n.is_finite() || residual_n.abs() > maximum_n {
        return Err(FiniteGapPlasticSensitivityError::BoundaryResidual {
            residual_n: residual_n.abs(),
            maximum_n,
        });
    }
    Ok(residual_n.abs())
}

fn pressure_force(
    model: &FiniteGapPlasticModel,
    pressure: &[f64],
    cx: &Cx<'_>,
) -> Result<f64, FiniteGapPlasticSensitivityError> {
    let mut sum = 0.0;
    let mut correction = 0.0;
    for chunk in pressure.chunks(64) {
        checkpoint(cx)?;
        for value in chunk {
            let corrected = value - correction;
            let next = sum + corrected;
            correction = (next - sum) - corrected;
            sum = next;
        }
    }
    let force = sum * (model.config.grid.cell_width_m * model.config.grid.cell_depth_m);
    if !force.is_finite() {
        return Err(FiniteGapPlasticSensitivityError::Invalid(
            "representable directional resultant force",
        ));
    }
    Ok(force)
}
