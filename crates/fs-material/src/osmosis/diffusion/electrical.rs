//! Floating electric fields for finite-rate, zero-current ionic diffusion.
//!
//! The original chemical force and one electric response use the SAME full
//! Maxwell--Stefan friction matrix. Eliminating the scalar field enforces
//! `sum z_i N_i=0` in the stationary-solvent frame; a final native solve checks
//! the combined force, its actual charge flux and its conjugate dissipation.
//! This is a local electroneutral, open-circuit limit. A network must separately
//! admit globally compatible potentials; arbitrary cyclic graphs cannot assume
//! zero current independently on every face. No Debye layer, Poisson field,
//! electrode reaction, imposed current, ionic activity model or heat is added.
//!
//! The dilute binary limit is the ambipolar/Nernst--Hartley law; its asymptotic
//! scope is discussed by Ghosal and Chen (2010),
//! <https://doi.org/10.1098/rspa.2010.0028>. The caller must qualify the original
//! virial chemical free energy and pair diffusivities for the assigned ions.
//! The force/friction, solvent-frame and conjugate current/entropy structure
//! follows Van-Brunt, Farrell and Monroe, equations 1.1--1.2, 6.5 and 16.13,
//! <https://arxiv.org/abs/2211.13505>. A common convex ionic face, shared with
//! source queries, preserves the neutral endpoint constraint. Individual
//! component logarithmic means would not preserve it for independent salts.

use super::{
    AmountConcentration, AqueousSolution, MAX_OSMOTIC_SOLUTES, MaxwellStefanWorkspaceV1, Result,
    SolutionDiffusionError, SolutionDiffusionReport, SolutionDiffusionSolve, SolutionMaxwellStefan,
    Sum, finite, logarithmic_mean_concentration, positive,
};
use fs_thermochem::{UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, stoichiometry::FARADAY_C_PER_MOL};

/// Explicit ionic identities and absolute physical charge budgets.
#[derive(Clone, Debug, PartialEq)]
pub struct ZeroCurrentDiffusion {
    valences: Box<[i32]>,
    maximum_charge_flux_mol_m2_s: f64,
    maximum_inventory_charge_mol: f64,
}

impl ZeroCurrentDiffusion {
    /// Original solute order; zero valence denotes a neutral co-solute. Both
    /// charge signs must be declared. No dissociation or missing ion is inferred.
    /// Charge units are moles of elementary-charge equivalents, not coulombs.
    pub fn new(
        valences: Vec<i32>,
        maximum_charge_flux_mol_m2_s: f64,
        maximum_inventory_charge_mol: f64,
    ) -> Result<Self> {
        if valences.len() > MAX_OSMOTIC_SOLUTES
            || !valences.iter().any(|&z| z > 0)
            || !valences.iter().any(|&z| z < 0)
            || !positive(maximum_charge_flux_mol_m2_s)
            || !positive(maximum_inventory_charge_mol)
        {
            return Err(SolutionDiffusionError::Invalid(
                "bounded opposite ion valences and positive absolute charge budgets required",
            ));
        }
        Ok(Self {
            valences: valences.into_boxed_slice(),
            maximum_charge_flux_mol_m2_s,
            maximum_inventory_charge_mol,
        })
    }
    pub fn valences(&self) -> &[i32] {
        &self.valences
    }
    pub fn heap_bytes(&self) -> usize {
        core::mem::size_of_val(self.valences.as_ref())
    }
    /// Check the actual finite inventory without charge repair or a relative
    /// allowance proportional to background salt. Returns its signed residual.
    pub fn admit_inventories(&self, amounts_mol: &[f64]) -> Result<f64> {
        if amounts_mol.len() != self.valences.len()
            || amounts_mol.iter().any(|n| !n.is_finite() || *n < 0.0)
        {
            return Err(SolutionDiffusionError::Invalid(
                "ionic inventory requires the original nonnegative solute layout",
            ));
        }
        let mut charge = Sum::default();
        for (&z, &amount) in self.valences.iter().zip(amounts_mol) {
            charge.add(finite(f64::from(z) * amount)?);
        }
        let residual_mol = charge.finish()?;
        if residual_mol.abs() > self.maximum_inventory_charge_mol {
            return Err(SolutionDiffusionError::InventoryCharge {
                residual_mol,
                maximum_mol: self.maximum_inventory_charge_mol,
            });
        }
        Ok(residual_mol)
    }

    /// A single convex interpolation preserves every linear charge constraint.
    /// Weight by reciprocal total ionic concentration along the straight path:
    /// its total is the logarithmic mean. Proportional salt profiles therefore
    /// retain their component logarithmic means exactly in real arithmetic.
    /// Neutral co-solutes keep their original logarithmic-mean face and exact
    /// zero-endpoint chain rule. One-sided-zero ions have singular chemical
    /// potentials under this reconstruction and explicitly refuse.
    pub(super) fn face_concentrations(&self, a: &[f64], b: &[f64], face: &mut [f64]) -> Result<()> {
        let mut totals = [Sum::default(), Sum::default()];
        for (i, &z) in self.valences.iter().enumerate() {
            if z != 0 {
                if (a[i] == 0.0) != (b[i] == 0.0) {
                    return Err(SolutionDiffusionError::SingularDiffusionPotential);
                }
                totals[0].add(a[i]);
                totals[1].add(b[i]);
            }
        }
        let [a_total, b_total] = totals.map(Sum::finish);
        let (a_total, b_total) = (a_total?, b_total?);
        let (lower, upper, low, high) = if a_total <= b_total {
            (a, b, a_total, b_total)
        } else {
            (b, a, b_total, a_total)
        };
        let theta = if low == high {
            0.5
        } else {
            let relative = (high - low) / low;
            if relative <= 1e-3 {
                // Stable divided difference of x/log(1+x) near equality.
                0.5 + relative
                    * (-1.0 / 12.0
                        + relative
                            * (1.0 / 24.0 + relative * (-19.0 / 720.0 + relative * 3.0 / 160.0)))
            } else {
                (logarithmic_mean_concentration(low, high)? - low) / (high - low)
            }
        };
        for (i, &z) in self.valences.iter().enumerate() {
            face[i] = if z == 0 {
                logarithmic_mean_concentration(a[i], b[i])?
            } else if lower[i] == upper[i] {
                lower[i]
            } else {
                finite((1.0 - theta) * lower[i] + theta * upper[i])?
            };
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn add_electric_force(
        &self,
        law: &SolutionMaxwellStefan,
        x: &[f64],
        drive: &mut [f64],
        concentration: f64,
        binary: &[f64],
        controls: SolutionDiffusionSolve,
        workspace: &mut MaxwellStefanWorkspaceV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<f64>> {
        let solvent = law.original_to_canonical[0];
        let current = |molar_flux: &[f64]| -> Result<f64> {
            let mut sum = Sum::default();
            for (i, &z) in self.valences.iter().enumerate() {
                let slot = law.original_to_canonical[i + 1];
                let flux = finite(molar_flux[slot] - (x[slot] / x[solvent]) * molar_flux[solvent])?;
                sum.add(finite(f64::from(z) * flux)?);
            }
            sum.finish().map_err(Into::into)
        };
        let chemical_current = current(
            law.owner
                .evaluate_molar_driving_forces(
                    x,
                    drive,
                    AmountConcentration::new(concentration),
                    binary,
                    controls.friction,
                    workspace,
                    &mut *cancelled,
                )?
                .molar_fluxes(),
        )?;
        // Normalize the electric response so very dilute carriers do not make
        // its driving vector vanish below the native solver's absolute scale.
        let mut electric = [0.0; MAX_OSMOTIC_SOLUTES + 1];
        let mut scale: f64 = 0.0;
        for (i, &z) in self.valences.iter().enumerate() {
            let slot = law.original_to_canonical[i + 1];
            electric[slot] = finite(x[slot] * f64::from(z))?;
            scale = scale.max(electric[slot].abs());
        }
        if scale == 0.0 {
            // A one-sided zero salt front with unequal mobilities may require
            // an infinite diffusion potential. Do not create trace carriers.
            return if chemical_current == 0.0 {
                Ok(None)
            } else {
                Err(SolutionDiffusionError::SingularDiffusionPotential)
            };
        }
        let mut balance = Sum::default();
        for &slot in &law.original_to_canonical[1..] {
            electric[slot] /= scale;
            balance.add(electric[slot]);
        }
        electric[solvent] = -balance.finish()?;
        let electric_current = current(
            law.owner
                .evaluate_molar_driving_forces(
                    x,
                    &electric[..x.len()],
                    AmountConcentration::new(concentration),
                    binary,
                    controls.friction,
                    workspace,
                    &mut *cancelled,
                )?
                .molar_fluxes(),
        )?;
        if electric_current >= 0.0 {
            return Err(SolutionDiffusionError::SingularDiffusionPotential);
        }
        let amplitude = finite(-chemical_current / electric_current)?;
        let gradient = finite(amplitude / scale)?;
        let mut balance = Sum::default();
        for &slot in &law.original_to_canonical[1..] {
            drive[slot] = finite(drive[slot] + finite(amplitude * electric[slot])?)?;
            balance.add(drive[slot]);
        }
        drive[solvent] = -balance.finish()?;
        Ok(Some(gradient))
    }

    pub(super) fn report(
        &self,
        law: &SolutionMaxwellStefan,
        molecular: SolutionDiffusionReport,
        gradient: Option<f64>,
        distance_m: f64,
        solvent_frame_fluxes: &[f64],
    ) -> Result<ZeroCurrentDiffusionReport> {
        let mut current = Sum::default();
        for (i, &z) in self.valences.iter().enumerate() {
            current.add(finite(
                f64::from(z) * solvent_frame_fluxes[law.original_to_canonical[i + 1]],
            )?);
        }
        let residual = current.finish()?;
        if residual.abs() > self.maximum_charge_flux_mol_m2_s {
            return Err(SolutionDiffusionError::ElectricCurrent {
                residual_mol_m2_s: residual,
                maximum_mol_m2_s: self.maximum_charge_flux_mol_m2_s,
            });
        }
        let difference = gradient.map(|g| finite(g * distance_m)).transpose()?;
        let voltage = difference
            .map(|difference| {
                finite(difference * (R * law.thermodynamics.temperature_k() / FARADAY_C_PER_MOL))
            })
            .transpose()?;
        Ok(ZeroCurrentDiffusionReport {
            molecular,
            dimensionless_potential_right_minus_left: difference,
            potential_right_minus_left_v: voltage,
            charge_flux_residual_mol_m2_s: residual,
            friction_evaluations: if gradient.is_some() { 3 } else { 2 },
        })
    }
}

/// Floating voltage uses no absolute potential datum. `None` means no finite
/// face carrier response determines it, with exactly zero uncorrected current.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZeroCurrentDiffusionReport {
    pub molecular: SolutionDiffusionReport,
    /// F (phi_right - phi_left) / (RT).
    pub dimensionless_potential_right_minus_left: Option<f64>,
    pub potential_right_minus_left_v: Option<f64>,
    pub charge_flux_residual_mol_m2_s: f64,
    /// Two response solves plus the final combined-force replay; two if no
    /// electric response exists. All use the original prepared workspace.
    pub friction_evaluations: usize,
}

impl SolutionMaxwellStefan {
    /// Original chemical/friction model with a floating zero-current field.
    /// Every finite endpoint must already meet the explicit charge budget.
    /// Successful evaluations allocate nothing; late refusals preserve output.
    #[allow(clippy::too_many_arguments)]
    pub fn flux_zero_current_into(
        &self,
        left: &AqueousSolution,
        right: &AqueousSolution,
        distance_m: f64,
        controls: SolutionDiffusionSolve,
        electrical: &ZeroCurrentDiffusion,
        workspace: &mut MaxwellStefanWorkspaceV1,
        output: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ZeroCurrentDiffusionReport> {
        self.flux_zero_current_with_binary_into(
            left,
            right,
            distance_m,
            controls,
            electrical,
            &self.binary_diffusivities_m2_s,
            workspace,
            output,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn flux_zero_current_with_binary_into(
        &self,
        left: &AqueousSolution,
        right: &AqueousSolution,
        distance_m: f64,
        controls: SolutionDiffusionSolve,
        electrical: &ZeroCurrentDiffusion,
        binary: &[f64],
        workspace: &mut MaxwellStefanWorkspaceV1,
        output: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ZeroCurrentDiffusionReport> {
        self.flux_with_binary_and_current_into(
            left,
            right,
            distance_m,
            controls,
            binary,
            Some(electrical),
            workspace,
            output,
            cancelled,
        )
        .map(|(_, electrical)| electrical.expect("ionic call retains its requested diagnostics"))
    }
}
