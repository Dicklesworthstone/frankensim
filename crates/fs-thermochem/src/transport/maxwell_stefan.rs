//! Neutral ideal-gas Maxwell–Stefan diffusion in explicitly chosen frames.
//!
//! With `z_i = j_i/M_i`, solve
//! `sum_j (x_j z_i - x_i z_j)/D_ij = -c [grad x_i + (x_i-Y_i) grad ln p]`
//! together with `sum_i M_i z_i = 0` in the mass-averaged frame or
//! `sum_i z_i = 0` in the molar-averaged frame. This is a coupled friction solve, not
//! independent Fick laws. It remains regular at exact zero mole fractions;
//! no division by `x_i` enters the linear system. Temperature specifies total
//! molar concentration, not a Soret force. Thermal diffusion, external
//! species forces, ions and nonideal activities are outside this closure.
//!
//! Primary formulation: COMSOL, *Multicomponent Gas Diffusion: Maxwell–Stefan
//! Description*, equations 3-75–3-80:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.chem/chem_ug_chemsptrans.08.257.html>.
//! This L1 module owns constitutive algebra and scratch storage, never evolution.

use crate::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_qty::{MolarMass, Pressure, Temperature};

/// Hard species/work bound for the dense constitutive solve.
pub const MAX_MAXWELL_STEFAN_SPECIES_V1: usize = 128;
/// Fixed operation-tree version of the barycentric friction solve.
pub const MAXWELL_STEFAN_EVALUATOR_VERSION_V1: u32 = 2;

/// Explicit checks on original constitutive equations and the mass frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaxwellStefanSolveV1 {
    /// Relative tolerance on each row's own driving/friction scale and net mass.
    pub relative_tolerance: f64,
    /// Absolute original-equation and gradient-compatibility tolerance [1/m].
    pub absolute_gradient_tolerance_per_m: f64,
    /// Absolute net diffusive mass-flux tolerance [kg/(m² s)].
    pub absolute_mass_flux_tolerance_kg_m2_s: f64,
}
impl Default for MaxwellStefanSolveV1 {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            absolute_gradient_tolerance_per_m: 1e-12,
            absolute_mass_flux_tolerance_kg_m2_s: 1e-16,
        }
    }
}

/// Original-equation and molar-frame checks for a molar-averaged evaluation.
/// A molar-frame flux generally transports nonzero net mass when species have
/// different molecular masses; it must not be used as a barycentric flux.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaxwellStefanMolarSolveV1 {
    /// Relative tolerance on each original row and the molar-frame sum.
    pub relative_tolerance: f64,
    /// Absolute original-equation and gradient-compatibility tolerance [1/m].
    pub absolute_gradient_tolerance_per_m: f64,
    /// Absolute sum of molar fluxes [mol/(m² s)].
    pub absolute_molar_flux_tolerance_mol_m2_s: f64,
}
impl Default for MaxwellStefanMolarSolveV1 {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            absolute_gradient_tolerance_per_m: 1e-12,
            absolute_molar_flux_tolerance_mol_m2_s: 1e-16,
        }
    }
}

#[derive(Clone, Copy)]
enum Frame {
    MassAveraged,
    MolarAveraged,
}

#[derive(Clone, Copy)]
struct Solve {
    relative_tolerance: f64,
    absolute_gradient_tolerance_per_m: f64,
    absolute_frame_tolerance: f64,
    frame: Frame,
}

/// Refusal before a complete borrowed evaluation can be published.
#[derive(Clone, Debug, PartialEq)]
pub enum MaxwellStefanErrorV1 {
    /// Invalid species, dimensions, fractions, coefficients or controls.
    Invalid(&'static str),
    /// A required algebraic intermediate is non-finite or unrepresentable.
    Arithmetic(&'static str),
    /// Scaled elimination found a zero or non-finite pivot.
    SingularSystem,
    /// An original equation or the mass-frame constraint missed its tolerance.
    Residual {
        /// Largest original Maxwell–Stefan row discrepancy [1/m].
        maximum_gradient_residual_per_m: f64,
        /// Signed sum of published species mass fluxes [kg/(m² s)].
        mass_frame_residual_kg_m2_s: f64,
    },
    /// An original equation or the molar-frame constraint missed its tolerance.
    MolarResidual {
        /// Largest original Maxwell–Stefan row discrepancy [1/m].
        maximum_gradient_residual_per_m: f64,
        /// Signed sum of published species molar fluxes [mol/(m² s)].
        molar_frame_residual_mol_m2_s: f64,
    },
    /// Bounded workspace reservation failed.
    AllocationRefused,
    /// Cancellation before complete publication.
    Cancelled,
}
impl core::fmt::Display for MaxwellStefanErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Maxwell-Stefan diffusion: {self:?}")
    }
}
impl std::error::Error for MaxwellStefanErrorV1 {}

/// Immutable canonical species and molecular masses for a diffusion closure.
/// Diffusivities remain explicit caller-supplied pair data at each actual state.
#[derive(Clone, Debug, PartialEq)]
pub struct MaxwellStefanV1 {
    species: Vec<SpeciesId>,
    molar_masses: Vec<MolarMass>,
    scaled_masses: Vec<f64>,
}
impl MaxwellStefanV1 {
    /// Canonicalize original identities and admit positive finite molar masses.
    /// Evaluation arrays follow the returned canonical `species()` order.
    pub fn new(mut species: Vec<(SpeciesId, MolarMass)>) -> Result<Self, MaxwellStefanErrorV1> {
        if !(1..=MAX_MAXWELL_STEFAN_SPECIES_V1).contains(&species.len())
            || species.iter().any(|(_, mass)| !positive(mass.value()))
        {
            return Err(MaxwellStefanErrorV1::Invalid("species count or molar mass"));
        }
        species.sort_by(|a, b| a.0.cmp(&b.0));
        if species.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(MaxwellStefanErrorV1::Invalid("duplicate species"));
        }
        let n = species.len();
        let maximum = species.iter().map(|(_, m)| m.value()).fold(0.0, f64::max);
        let mut ids = reserved(n)?;
        let mut masses = reserved(n)?;
        let mut scaled = reserved(n)?;
        for (id, mass) in species {
            let value = mass.value() / maximum;
            if !positive(value) {
                return Err(MaxwellStefanErrorV1::Arithmetic("molar-mass scaling"));
            }
            ids.push(id);
            masses.push(mass);
            scaled.push(value);
        }
        Ok(Self {
            species: ids,
            molar_masses: masses,
            scaled_masses: scaled,
        })
    }
    /// Original species identities in deterministic canonical order.
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }
    /// Original molar masses [kg/mol] in canonical species order.
    pub fn molar_masses(&self) -> &[MolarMass] {
        &self.molar_masses
    }
    /// Allocate the complete reusable matrix, solve and output workspace once.
    pub fn workspace(&self) -> Result<MaxwellStefanWorkspaceV1, MaxwellStefanErrorV1> {
        let n = self.species.len();
        Ok(MaxwellStefanWorkspaceV1 {
            matrix: filled(n * n, 0.0)?,
            ratios: filled(n * n, 0.0)?,
            right: filled(n, 0.0)?,
            row_scales: filled(n, 0.0)?,
            driving: filled(n, 0.0)?,
            molar_flux: filled(n, 0.0)?,
            mass_flux: filled(n, 0.0)?,
        })
    }

    /// Evaluate one spatial direction without successful-path allocation.
    ///
    /// `binary_diffusivities_m2_s` is an exactly symmetric row-major `n*n`
    /// matrix with positive finite off-diagonals and zero diagonal. Its source
    /// temperature, pressure and validity are the caller's responsibility.
    /// Fractions may contain exact zeros; neither fractions nor fluxes are
    /// clipped or projected. Every original friction row, including the one
    /// replaced by the barycentric constraint, is checked before return.
    /// Naturally rounded subnormal or zero fluxes must pass those same checks;
    /// a small driving force alone is never grounds to suppress a response.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate<'a>(
        &self,
        mole_fractions: &[f64],
        gradients_per_m: &[f64],
        temperature: Temperature,
        pressure: Pressure,
        log_pressure_gradient_per_m: f64,
        binary_diffusivities_m2_s: &[f64],
        controls: MaxwellStefanSolveV1,
        workspace: &'a mut MaxwellStefanWorkspaceV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanEvaluationV1<'a>, MaxwellStefanErrorV1> {
        self.evaluate_in_frame(
            mole_fractions,
            gradients_per_m,
            temperature,
            pressure,
            log_pressure_gradient_per_m,
            binary_diffusivities_m2_s,
            Solve {
                relative_tolerance: controls.relative_tolerance,
                absolute_gradient_tolerance_per_m: controls.absolute_gradient_tolerance_per_m,
                absolute_frame_tolerance: controls.absolute_mass_flux_tolerance_kg_m2_s,
                frame: Frame::MassAveraged,
            },
            workspace,
            cancelled,
        )
    }

    /// Evaluate in the molar-averaged frame with zero summed molar flux.
    ///
    /// The original friction equations, canonical source order and pressure
    /// driving force are identical to [`Self::evaluate`]. Only the reference
    /// velocity changes: `sum(z_i)=0` replaces `sum(M_i z_i)=0` in the same
    /// scaled solve. Molar flux is the primary rounded output in this frame;
    /// no conversion through a potentially underflowed mass flux is required.
    /// Exact-zero fractions, workspace reuse and cancellation have the same
    /// contract. Thermal diffusion and nonideal activities remain excluded.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_molar<'a>(
        &self,
        mole_fractions: &[f64],
        gradients_per_m: &[f64],
        temperature: Temperature,
        pressure: Pressure,
        log_pressure_gradient_per_m: f64,
        binary_diffusivities_m2_s: &[f64],
        controls: MaxwellStefanMolarSolveV1,
        workspace: &'a mut MaxwellStefanWorkspaceV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanMolarEvaluationV1<'a>, MaxwellStefanErrorV1> {
        self.evaluate_in_frame(
            mole_fractions,
            gradients_per_m,
            temperature,
            pressure,
            log_pressure_gradient_per_m,
            binary_diffusivities_m2_s,
            Solve {
                relative_tolerance: controls.relative_tolerance,
                absolute_gradient_tolerance_per_m: controls.absolute_gradient_tolerance_per_m,
                absolute_frame_tolerance: controls.absolute_molar_flux_tolerance_mol_m2_s,
                frame: Frame::MolarAveraged,
            },
            workspace,
            cancelled,
        )
        .map(|evaluation| MaxwellStefanMolarEvaluationV1 { evaluation })
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_in_frame<'a>(
        &self,
        mole_fractions: &[f64],
        gradients_per_m: &[f64],
        temperature: Temperature,
        pressure: Pressure,
        log_pressure_gradient_per_m: f64,
        binary_diffusivities_m2_s: &[f64],
        controls: Solve,
        workspace: &'a mut MaxwellStefanWorkspaceV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanEvaluationV1<'a>, MaxwellStefanErrorV1> {
        poll(&mut cancelled)?;
        let n = self.species.len();
        if mole_fractions.len() != n
            || gradients_per_m.len() != n
            || binary_diffusivities_m2_s.len() != n * n
            || workspace.right.len() != n
            || !positive(temperature.value())
            || !positive(pressure.value())
            || !log_pressure_gradient_per_m.is_finite()
            || !positive(controls.relative_tolerance)
            || controls.relative_tolerance > 1e-3
            || !positive(controls.absolute_gradient_tolerance_per_m)
            || !positive(controls.absolute_frame_tolerance)
            || mole_fractions
                .iter()
                .any(|x| !x.is_finite() || !(0.0..=1.0).contains(x))
            || gradients_per_m.iter().any(|x| !x.is_finite())
        {
            return Err(MaxwellStefanErrorV1::Invalid(
                "state, array dimensions or numerical controls",
            ));
        }
        let fraction_sum = sum(mole_fractions.iter().copied());
        if (fraction_sum - 1.0).abs() > 64.0 * f64::EPSILON * n as f64 {
            return Err(MaxwellStefanErrorV1::Invalid(
                "mole fractions do not close the simplex",
            ));
        }
        let gradient_sum = sum(gradients_per_m.iter().copied());
        let gradient_scale = sum(gradients_per_m.iter().map(|g| g.abs()));
        if !gradient_scale.is_finite()
            || gradient_sum.abs()
                > controls.absolute_gradient_tolerance_per_m
                    + controls.relative_tolerance * gradient_scale
        {
            return Err(MaxwellStefanErrorV1::Invalid(
                "mole-fraction gradients do not sum to zero",
            ));
        }
        let mean_mass = sum(mole_fractions
            .iter()
            .zip(&self.molar_masses)
            .map(|(x, m)| x * m.value()));
        let concentration = (pressure.value() / R) / temperature.value();
        if !positive(mean_mass) || !positive(concentration) {
            return Err(MaxwellStefanErrorV1::Arithmetic(
                "mean molar mass or total concentration",
            ));
        }
        let mut minimum_diffusion = f64::INFINITY;
        for i in 0..n {
            poll(&mut cancelled)?;
            if binary_diffusivities_m2_s[i * n + i] != 0.0 {
                return Err(MaxwellStefanErrorV1::Invalid(
                    "nonzero binary-diffusion diagonal",
                ));
            }
            for j in i + 1..n {
                let value = binary_diffusivities_m2_s[i * n + j];
                if !positive(value) || value != binary_diffusivities_m2_s[j * n + i] {
                    return Err(MaxwellStefanErrorV1::Invalid(
                        "binary diffusion is not positive symmetric",
                    ));
                }
                minimum_diffusion = minimum_diffusion.min(value);
            }
        }
        let mut drive_scale = 0.0_f64;
        let mut gauge = 0;
        for i in 0..n {
            let mass_fraction = (mole_fractions[i] * self.molar_masses[i].value()) / mean_mass;
            let drive = gradients_per_m[i]
                + (mole_fractions[i] - mass_fraction) * log_pressure_gradient_per_m;
            if !drive.is_finite() {
                return Err(MaxwellStefanErrorV1::Arithmetic(
                    "pressure/composition driving force",
                ));
            }
            workspace.driving[i] = drive;
            drive_scale = drive_scale.max(drive.abs());
            if mole_fractions[i] > mole_fractions[gauge] {
                gauge = i;
            }
        }
        workspace.mass_flux.fill(0.0);
        workspace.molar_flux.fill(0.0);
        if n == 1 || drive_scale == 0.0 {
            if drive_scale > controls.absolute_gradient_tolerance_per_m {
                return Err(MaxwellStefanErrorV1::Invalid(
                    "single-species nonzero diffusion force",
                ));
            }
            poll(&mut cancelled)?;
            return Ok(MaxwellStefanEvaluationV1 {
                mass_fluxes: &workspace.mass_flux,
                molar_fluxes: &workspace.molar_flux,
                maximum_gradient_residual_per_m: drive_scale,
                frame_residual: 0.0,
                drag_entropy_production_w_m3_k: Some(0.0),
            });
        }
        workspace.matrix.fill(0.0);
        workspace.ratios.fill(0.0);
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let ratio = minimum_diffusion / binary_diffusivities_m2_s[i * n + j];
                if !positive(ratio) {
                    return Err(MaxwellStefanErrorV1::Arithmetic("binary-diffusivity ratio"));
                }
                workspace.ratios[i * n + j] = ratio;
                workspace.matrix[i * n + i] += mole_fractions[j] * ratio;
                workspace.matrix[i * n + j] = -mole_fractions[i] * ratio;
            }
            workspace.right[i] = -workspace.driving[i] / drive_scale;
        }
        match controls.frame {
            Frame::MassAveraged => {
                workspace.matrix[gauge * n..(gauge + 1) * n].copy_from_slice(&self.scaled_masses)
            }
            Frame::MolarAveraged => workspace.matrix[gauge * n..(gauge + 1) * n].fill(1.0),
        }
        workspace.right[gauge] = 0.0;
        solve(
            &mut workspace.matrix,
            &mut workspace.right,
            &mut workspace.row_scales,
            n,
            &mut cancelled,
        )?;
        for i in 0..n {
            // The common c*D_min*drive scale can underflow even when the
            // solved amplitude makes the actual species flux representable.
            // Include every factor before rounding into the binary64 range.
            let molecular_mass = self.molar_masses[i].value();
            let (mass, molar) = match controls.frame {
                Frame::MassAveraged => {
                    let mass = scaled_product_ratio(
                        &[
                            concentration,
                            minimum_diffusion,
                            drive_scale,
                            workspace.right[i],
                            molecular_mass,
                        ],
                        &[],
                    )?;
                    (mass, mass / molecular_mass)
                }
                Frame::MolarAveraged => {
                    let molar = scaled_product_ratio(
                        &[
                            concentration,
                            minimum_diffusion,
                            drive_scale,
                            workspace.right[i],
                        ],
                        &[],
                    )?;
                    (molar * molecular_mass, molar)
                }
            };
            if !molar.is_finite() || !mass.is_finite() {
                return Err(MaxwellStefanErrorV1::Arithmetic("species diffusion flux"));
            }
            workspace.mass_flux[i] = mass;
            workspace.molar_flux[i] = molar;
            // Residuals concern the primary flux actually published in the
            // chosen frame, including conversion rounding in mass mode.
            workspace.right[i] = match controls.frame {
                Frame::MassAveraged => scaled_product_ratio(
                    &[mass],
                    &[
                        molecular_mass,
                        concentration,
                        minimum_diffusion,
                        drive_scale,
                    ],
                )?,
                Frame::MolarAveraged => scaled_product_ratio(
                    &[molar],
                    &[concentration, minimum_diffusion, drive_scale],
                )?,
            };
        }
        let constrained_flux = match controls.frame {
            Frame::MassAveraged => &workspace.mass_flux,
            Frame::MolarAveraged => &workspace.molar_flux,
        };
        let frame_residual = sum(constrained_flux.iter().copied());
        let frame_scale = sum(constrained_flux.iter().map(|j| j.abs()));
        let mut maximum_residual = 0.0_f64;
        let mut rows_close = true;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut terms = [0.0; MAX_MAXWELL_STEFAN_SPECIES_V1];
            let wi = workspace.right[i];
            for j in 0..n {
                if i != j {
                    let wj = workspace.right[j];
                    terms[j] = (mole_fractions[j] * wi - mole_fractions[i] * wj)
                        * workspace.ratios[i * n + j]
                        * drive_scale;
                }
            }
            let residual = workspace.driving[i] + sum(terms[..n].iter().copied());
            let row_scale = workspace.driving[i].abs() + sum(terms[..n].iter().map(|t| t.abs()));
            if !residual.is_finite() || !row_scale.is_finite() {
                return Err(MaxwellStefanErrorV1::Arithmetic(
                    "original diffusion residual",
                ));
            }
            maximum_residual = maximum_residual.max(residual.abs());
            rows_close &= residual.abs()
                <= controls.absolute_gradient_tolerance_per_m
                    + controls.relative_tolerance * row_scale;
        }
        if !frame_residual.is_finite() || !frame_scale.is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic("diffusion-frame residual"));
        }
        if !rows_close
            || frame_residual.abs()
                > controls.absolute_frame_tolerance + controls.relative_tolerance * frame_scale
        {
            return Err(match controls.frame {
                Frame::MassAveraged => MaxwellStefanErrorV1::Residual {
                    maximum_gradient_residual_per_m: maximum_residual,
                    mass_frame_residual_kg_m2_s: frame_residual,
                },
                Frame::MolarAveraged => MaxwellStefanErrorV1::MolarResidual {
                    maximum_gradient_residual_per_m: maximum_residual,
                    molar_frame_residual_mol_m2_s: frame_residual,
                },
            });
        }
        let entropy = drag_entropy(
            mole_fractions,
            &workspace.molar_flux,
            binary_diffusivities_m2_s,
            concentration,
        );
        poll(&mut cancelled)?;
        Ok(MaxwellStefanEvaluationV1 {
            mass_fluxes: &workspace.mass_flux,
            molar_fluxes: &workspace.molar_flux,
            maximum_gradient_residual_per_m: maximum_residual,
            frame_residual,
            drag_entropy_production_w_m3_k: entropy,
        })
    }
}

/// Reusable bounded constitutive scratch, with no evolving physical state.
#[derive(Clone, Debug)]
pub struct MaxwellStefanWorkspaceV1 {
    matrix: Vec<f64>,
    ratios: Vec<f64>,
    right: Vec<f64>,
    row_scales: Vec<f64>,
    driving: Vec<f64>,
    molar_flux: Vec<f64>,
    mass_flux: Vec<f64>,
}

impl MaxwellStefanWorkspaceV1 {
    /// Allocated numerical array capacities in bytes. Immutable species and
    /// fixed-size owner values are excluded. Higher-level spatial operators
    /// include this original solve storage in their prepared memory budgets.
    pub fn workspace_bytes(&self) -> usize {
        [
            &self.matrix,
            &self.ratios,
            &self.right,
            &self.row_scales,
            &self.driving,
            &self.molar_flux,
            &self.mass_flux,
        ]
        .iter()
        .map(|values| values.capacity() * core::mem::size_of::<f64>())
        .sum()
    }
}

/// Complete borrowed diffusive fluxes and actual equation residuals.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellStefanEvaluationV1<'a> {
    mass_fluxes: &'a [f64],
    molar_fluxes: &'a [f64],
    maximum_gradient_residual_per_m: f64,
    frame_residual: f64,
    drag_entropy_production_w_m3_k: Option<f64>,
}
impl MaxwellStefanEvaluationV1<'_> {
    /// Canonical diffusive mass flux relative to the mass-averaged velocity.
    pub const fn mass_fluxes(&self) -> &[f64] {
        self.mass_fluxes
    }
    /// Canonical molar flux in that same mass-averaged frame [mol/(m² s)].
    /// Its sum is generally nonzero when molecular masses differ.
    pub const fn molar_fluxes(&self) -> &[f64] {
        self.molar_fluxes
    }
    /// Largest absolute residual among all original friction equations [1/m].
    pub const fn maximum_gradient_residual_per_m(&self) -> f64 {
        self.maximum_gradient_residual_per_m
    }
    /// Signed sum of the published mass fluxes [kg/(m² s)].
    pub const fn mass_frame_residual_kg_m2_s(&self) -> f64 {
        self.frame_residual
    }
    /// Local nonnegative pair-friction entropy production [W/(m³ K)].
    /// None denotes a nonzero flux at an exact-zero composition boundary or an
    /// unrepresentable diagnostic. This is not a discrete transport inequality.
    pub const fn drag_entropy_production_w_m3_k(&self) -> Option<f64> {
        self.drag_entropy_production_w_m3_k
    }
}

/// Complete fluxes relative to the molar-averaged velocity. This distinct type
/// cannot silently substitute for a mass-averaged constitutive evaluation.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellStefanMolarEvaluationV1<'a> {
    evaluation: MaxwellStefanEvaluationV1<'a>,
}
impl MaxwellStefanMolarEvaluationV1<'_> {
    /// Canonical molar fluxes [mol/(m² s)] whose sum closes the declared budget.
    pub const fn molar_fluxes(&self) -> &[f64] {
        self.evaluation.molar_fluxes
    }
    /// Canonical mass fluxes in the same molar frame [kg/(m² s)]. Their sum is
    /// generally nonzero; this is not the barycentric momentum-frame flux.
    pub const fn mass_fluxes(&self) -> &[f64] {
        self.evaluation.mass_fluxes
    }
    /// Largest original Maxwell–Stefan row discrepancy [1/m].
    pub const fn maximum_gradient_residual_per_m(&self) -> f64 {
        self.evaluation.maximum_gradient_residual_per_m
    }
    /// Signed sum of published species molar fluxes [mol/(m² s)].
    pub const fn molar_frame_residual_mol_m2_s(&self) -> f64 {
        self.evaluation.frame_residual
    }
    /// Frame-invariant pair-friction entropy production [W/(m³ K)].
    /// The same exact-zero-composition and representability boundaries apply
    /// as to the mass-averaged evaluation; no discrete entropy proof is implied.
    pub const fn drag_entropy_production_w_m3_k(&self) -> Option<f64> {
        self.evaluation.drag_entropy_production_w_m3_k
    }
}

fn solve(
    a: &mut [f64],
    b: &mut [f64],
    scales: &mut [f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MaxwellStefanErrorV1> {
    for i in 0..n {
        scales[i] = a[i * n..(i + 1) * n]
            .iter()
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        if !positive(scales[i]) {
            return Err(MaxwellStefanErrorV1::SingularSystem);
        }
    }
    for k in 0..n {
        poll(cancelled)?;
        let mut pivot = k;
        for i in k + 1..n {
            if a[i * n + k].abs() / scales[i] > a[pivot * n + k].abs() / scales[pivot] {
                pivot = i;
            }
        }
        if a[pivot * n + k] == 0.0 || !a[pivot * n + k].is_finite() {
            return Err(MaxwellStefanErrorV1::SingularSystem);
        }
        if pivot != k {
            for j in 0..n {
                a.swap(k * n + j, pivot * n + j);
            }
            b.swap(k, pivot);
            scales.swap(k, pivot);
        }
        for i in k + 1..n {
            let factor = a[i * n + k] / a[k * n + k];
            a[i * n + k] = 0.0;
            for j in k + 1..n {
                a[i * n + j] -= factor * a[k * n + j];
            }
            b[i] -= factor * b[k];
        }
    }
    for i in (0..n).rev() {
        poll(cancelled)?;
        let rhs = b[i] - sum((i + 1..n).map(|j| a[i * n + j] * b[j]));
        b[i] = rhs / a[i * n + i];
        if !b[i].is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic("diffusion linear solve"));
        }
    }
    Ok(())
}

fn drag_entropy(x: &[f64], z: &[f64], diffusion: &[f64], concentration: f64) -> Option<f64> {
    let n = x.len();
    if x.iter().zip(z).any(|(x, z)| *x == 0.0 && *z != 0.0) {
        return None;
    }
    let mut total = 0.0;
    for i in 0..n {
        if x[i] == 0.0 {
            continue;
        }
        for j in i + 1..n {
            if x[j] == 0.0 {
                continue;
            }
            let velocity = (z[i] / concentration) / x[i] - (z[j] / concentration) / x[j];
            let term =
                (R * concentration * x[i] * x[j] / diffusion[i * n + j]) * velocity * velocity;
            if !term.is_finite() {
                return None;
            }
            total += term;
        }
    }
    total.is_finite().then_some(total)
}
/// Balance binary exponents before final rounding, without forming a potentially
/// unrepresentable common scale. A rounded zero is only a candidate: the caller
/// still checks the original constitutive equations with its unchanged budgets.
fn scaled_product_ratio(
    numerator: &[f64],
    denominator: &[f64],
) -> Result<f64, MaxwellStefanErrorV1> {
    let mut direct = 1.0;
    let mut ordinary = true;
    let mut negative = false;
    let mut zero = false;
    for &factor in numerator {
        if !factor.is_finite() {
            return Err(MaxwellStefanErrorV1::Arithmetic("scaled diffusion product"));
        }
        negative ^= factor.is_sign_negative();
        zero |= factor == 0.0;
        direct *= factor;
        ordinary &= direct.is_normal();
    }
    for &factor in denominator {
        if !positive(factor) {
            return Err(MaxwellStefanErrorV1::Arithmetic(
                "scaled diffusion quotient",
            ));
        }
        direct /= factor;
        ordinary &= direct.is_normal();
    }
    if zero {
        return Ok(if negative { -0.0 } else { 0.0 });
    }
    if ordinary {
        return Ok(direct);
    }
    let mut mantissa = 1.0;
    let mut exponent = 0_i32;
    for &factor in numerator {
        let (part, power) = binary_parts(factor.abs());
        mantissa *= part;
        exponent += power;
        if mantissa >= 2.0 {
            mantissa *= 0.5;
            exponent += 1;
        }
    }
    for &factor in denominator {
        let (part, power) = binary_parts(factor);
        mantissa /= part;
        exponent -= power;
        if mantissa < 1.0 {
            mantissa *= 2.0;
            exponent -= 1;
        }
    }
    let magnitude = if exponent > 1023 {
        return Err(MaxwellStefanErrorV1::Arithmetic(
            "scaled diffusion result overflow",
        ));
    } else if exponent >= -1022 {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent >= -1075 {
        // Both factors are normal until the last multiplication, so gradual
        // underflow is rounded once (including the half-minimum ties-to-even).
        (mantissa * f64::from_bits(((exponent + 2045) as u64) << 52)) * f64::MIN_POSITIVE
    } else {
        0.0
    };
    if !magnitude.is_finite() {
        return Err(MaxwellStefanErrorV1::Arithmetic(
            "scaled diffusion result overflow",
        ));
    }
    Ok(if negative { -magnitude } else { magnitude })
}

/// Exact decomposition of positive normal or subnormal binary64 into m*2^e,
/// with 1 <= m < 2. No logarithm or reciprocal of a tiny input is needed.
fn binary_parts(value: f64) -> (f64, i32) {
    let bits = value.to_bits();
    let power = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & 0x000f_ffff_ffff_ffff;
    if power != 0 {
        (f64::from_bits((1023_u64 << 52) | fraction), power - 1023)
    } else {
        let highest = 63 - fraction.leading_zeros();
        (
            fraction as f64 / (1_u64 << highest) as f64,
            highest as i32 - 1074,
        )
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaxwellStefanErrorV1> {
    if cancelled() {
        Err(MaxwellStefanErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn reserved<T>(count: usize) -> Result<Vec<T>, MaxwellStefanErrorV1> {
    let mut v = Vec::new();
    v.try_reserve_exact(count)
        .map_err(|_| MaxwellStefanErrorV1::AllocationRefused)?;
    Ok(v)
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, MaxwellStefanErrorV1> {
    let mut v = reserved(count)?;
    v.resize(count, value);
    Ok(v)
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        let y = value - correction;
        let next = total + y;
        correction = (next - total) - y;
        total = next;
    }
    total
}

#[cfg(test)]
mod tests;
