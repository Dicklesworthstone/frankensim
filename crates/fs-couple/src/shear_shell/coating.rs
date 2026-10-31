//! Perfectly bonded films on the original material shell reference surface.
//!
//! Storage is frozen at an authored frequency with an enforced drift bound.
//! Each film's measured loss is read at every retained mode's own frequency;
//! a positive diagonal congruence preserves modal cross-loss passivity. This
//! geometric-mean off-diagonal reading is a modal approximation, not a full
//! hereditary time-domain shell or a frequency-dependent eigenproblem.
//! The explicitly selected [`relaxation`] path instead builds equilibrium
//! sections and projects the complete passive coating memory on that basis.

use super::{
    AcousticRealizeError, PlateMaterialModel, ResolvedShearShell, compile_material_shear_shell,
    invalid, owner, plate_property, poll,
};
use fs_exec::CancelGate;
use fs_material::coating::{CoatingFrequencySlice, ResolvedCoating};
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use fs_material::sheet::SheetMatrix;
use fs_plate::shell::mitc::{
    MitcBudget, NodalDirectors, ShearShellSection,
    viscosity::{ShearShellViscosity, project_viscosity, projection_terms},
};
use fs_qty::Pressure;

/// Equilibrium sections and source-fitted passive Maxwell memory.
pub mod relaxation;
/// In-band source-law structural response and uncoated body comparison.
pub mod harmonic;
pub mod strain;

/// Cold-work ceiling for the complete coating stack across one shell.
pub const MAX_SHELL_COATINGS: usize = 64;
/// Bound retained per-film, per-facet coefficient records before allocation.
pub const MAX_COATING_SECTION_ENTRIES: usize = 262_144;

/// Side of the substrate with respect to its oriented directors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoatingSide {
    /// Positive director side, starting at +substrate_thickness/2.
    Positive,
    /// Negative director side, starting at -substrate_thickness/2.
    Negative,
}

/// One completely specified, source-backed surface layer.
#[derive(Clone, Debug)]
pub struct ShellCoating {
    /// Existing material region to coat completely; no triangle overlap is inferred.
    pub region: String,
    /// Original measured spectrum, source state and authored isotropic lift.
    /// Its strain limit remains a small-strain modeling assumption; the linear
    /// modal runtime does not enforce local displacement-amplitude bounds.
    pub material: ResolvedCoating,
    /// Physical side, independent of exterior/cavity pressure ownership.
    pub side: CoatingSide,
    /// Actual added film thickness [m]; substrate thickness is unchanged.
    pub thickness_m: f64,
    /// Storage reference [rad/s], explicitly inside both material source bands.
    pub reference_omega_rad_s: f64,
    /// Maximum |E'(reference)-E'(mode)|/E'(mode), finite in [0,1).
    /// Every retained mode must pass this authored frozen-storage allowance.
    pub max_storage_relative_drift: f64,
}

/// Original film input and its exact placement about the substrate midplane.
#[derive(Clone, Debug)]
pub struct ResolvedShellCoating {
    input: ShellCoating,
    area_m2: f64,
    mass_kg: f64,
    omega_band: (f64, f64),
    z_bottom_m: f64,
    z_top_m: f64,
    reference_storage_pa: f64,
    reference_viscosity_pa_s: f64,
}

impl ResolvedShellCoating {
    /// Unchanged material source, geometry and approximation policy.
    pub const fn input(&self) -> &ShellCoating {
        &self.input
    }
    /// Coated reference-surface area [m²].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Added physical film mass [kg], excluding substrate and other layers.
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Intersected substrate/film support [rad/s]; storage drift is checked separately.
    pub const fn omega_band_rad_s(&self) -> (f64, f64) {
        self.omega_band
    }
    /// Lower and upper film faces [m], measured from the substrate midplane.
    pub const fn faces_m(&self) -> (f64, f64) {
        (self.z_bottom_m, self.z_top_m)
    }
    /// Actual basis storage [Pa]: source E'(reference) for a slice, or the
    /// declared model's E0 when hereditary coating dynamics are selected.
    pub const fn reference_storage_pa(&self) -> f64 {
        self.reference_storage_pa
    }

    fn loss_scale(&self, omega: f64) -> Result<f64, AcousticRealizeError> {
        if !omega.is_finite()
            || omega <= 0.0
            || omega < self.omega_band.0
            || omega > self.omega_band.1
        {
            return Err(invalid(
                "retained mode lies outside the coating/substrate source support",
            ));
        }
        let (storage, loss) = self.input.material.law().modulus(omega).map_err(owner)?;
        let drift = (storage - self.reference_storage_pa).abs() / storage;
        if !drift.is_finite() || drift > self.input.max_storage_relative_drift {
            return Err(owner(format_args!(
                "coating on region {} exceeds its frozen-storage allowance at {} rad/s: relative drift {}, allowed {}",
                self.input.region, omega, drift, self.input.max_storage_relative_drift
            )));
        }
        let viscosity = loss / omega;
        if self.reference_viscosity_pa_s == 0.0 {
            return if viscosity == 0.0 {
                Ok(1.0)
            } else {
                Err(invalid(
                    "zero reference film loss cannot represent a nonzero modal loss",
                ))
            };
        }
        let scale = fs_math::det::sqrt(viscosity / self.reference_viscosity_pa_s);
        if !scale.is_finite() {
            return Err(invalid(
                "source coating modal loss ratio is unrepresentable",
            ));
        }
        Ok(scale)
    }
}

impl ResolvedShearShell {
    /// Original coating records in caller order. Within each region and side,
    /// inputs are ordered from the substrate outward, including ground layers.
    pub fn coatings(&self) -> &[ResolvedShellCoating] {
        &self.coatings
    }

    /// Replace the complete coating stack transactionally. Empty removes it.
    /// Regions retain their original substrate mass/thickness policies; films
    /// add mass and do not consume a fixed-mass substrate allowance. Geometry,
    /// material axes, pressure traces and authored attachment offsets survive.
    /// Perfect bonding is assumed; curing stress, penetration and slip are absent.
    pub fn with_coatings(
        &self,
        inputs: Vec<ShellCoating>,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        admit_extent(inputs.len(), self.mesh.tris.len())?;
        let substrate = compile_material_shear_shell(
            self.mesh.clone(),
            self.regions.iter().map(|r| r.input.clone()).collect(),
            self.budget,
            gate,
        )?;
        let coated = apply_coatings(substrate, inputs, self.coating_relaxation, gate)?;
        let mut coated = coated.with_interfaces(
            self.interfaces.iter().map(|i| i.input().clone()).collect(),
            self.interface_budget,
            gate,
        )?;
        coated.pressure_surfaces = self.pressure_surfaces.clone();
        poll(gate)?;
        Ok(coated)
    }

    /// Check source support and the authored frozen-storage bound at actual
    /// physical frequencies. Success does not assert a hereditary shell model.
    pub fn validate_coating_frequencies(
        &self,
        omega: &[f64],
        gate: &CancelGate,
    ) -> Result<(), AcousticRealizeError> {
        for film in &self.coatings {
            for &w in omega {
                poll(gate)?;
                if self.coating_relaxation.is_some() {
                    if !w.is_finite() || w < film.omega_band.0 || w > film.omega_band.1 {
                        return Err(invalid("retained hereditary shell mode lies outside coating/source fit support"));
                    }
                } else {
                    film.loss_scale(w)?;
                }
            }
        }
        poll(gate)
    }

    pub(super) fn loss_projection_terms(&self, modes: usize) -> Option<usize> {
        let distinct: std::collections::BTreeSet<_> = self
            .regions
            .iter()
            .map(|r| r.loss_reference_omega.map(f64::to_bits))
            .collect();
        let projected = projection_terms(self.mesh.nodes.len(), self.mesh.tris.len(), modes)?
            .checked_mul(distinct.len().checked_add(self.coatings.len())?)?;
        // In the hereditary path the film forms are projected once instead of
        // frozen film C. Include complete branch accumulation and factor work.
        let memory = self.coating_relaxation.map_or(Some(0), |options|
            options.projection_terms(modes, self.coatings.len()))?;
        projected.checked_add(memory)
    }

    /// Project all substrate and film loss on the SAME physical mode vectors.
    /// Anchored wood loss uses its original law. Each film instead uses its own
    /// E''(omega)/omega, exactly on modal diagonals and by geometric-mean
    /// congruence on cross terms. Storage and condensation remain at the authored
    /// reference; no extra damping, modal diagonalization or memory is inserted.
    /// Returns row-major modal C [1/s] when vectors are M-orthonormal.
    #[allow(clippy::too_many_arguments)]
    pub fn project_material_viscosity(
        &self,
        directors: &NodalDirectors,
        physical_modes: &[Vec<f64>],
        omega: &[f64],
        budget: MitcBudget,
        max_projection_terms: usize,
        gate: &CancelGate,
    ) -> Result<Vec<f64>, AcousticRealizeError> {
        poll(gate)?;
        let n = physical_modes.len();
        if n == 0
            || n > super::MAX_SHELL_MODES
            || omega.len() != n
            || omega.iter().any(|w| {
                !w.is_finite()
                    || *w <= 0.0
                    || self
                        .regions
                        .iter()
                        .any(|r| *w < r.omega_band.0 || *w > r.omega_band.1)
            })
            || self
                .loss_projection_terms(n)
                .is_none_or(|w| w > max_projection_terms)
        {
            return Err(invalid(
                "material shell loss needs supported modes and complete per-law projection allowance",
            ));
        }
        self.validate_coating_frequencies(omega, gate)?;
        let spectrum = super::constant_loss::ConstantLossSpectrum::fit();
        let mut matrix = Vec::new();
        for (law, viscosity) in self.loss_groups()? {
            poll(gate)?;
            let mut part = project_viscosity(
                &self.mesh,
                directors,
                &self.sections,
                &viscosity,
                physical_modes,
                budget,
                max_projection_terms,
                || gate.is_requested(),
            )
            .map_err(owner)?;
            if let Some(reference) = law {
                super::constant_loss::rescale(&mut part, omega, reference, &spectrum);
            }
            if matrix.is_empty() {
                matrix = part;
            } else {
                for (out, value) in matrix.iter_mut().zip(part) {
                    *out += value;
                }
            }
        }
        for (film, viscosity) in self.coatings.iter().zip(&self.coating_viscosity)
            .filter(|_| self.coating_relaxation.is_none()) {
            poll(gate)?;
            let scales: Vec<_> = omega
                .iter()
                .map(|&w| film.loss_scale(w))
                .collect::<Result<_, _>>()?;
            let part = project_viscosity(
                &self.mesh,
                directors,
                &self.sections,
                viscosity,
                physical_modes,
                budget,
                max_projection_terms,
                || gate.is_requested(),
            )
            .map_err(owner)?;
            for i in 0..n {
                for j in i..n {
                    let value = part[n * i + j] * (scales[i] * scales[j]);
                    matrix[n * i + j] += value;
                    if i != j {
                        matrix[n * j + i] += value;
                    }
                }
            }
        }
        if matrix.iter().any(|v| !v.is_finite()) {
            return Err(invalid("combined shell material loss is unrepresentable"));
        }
        poll(gate)?;
        Ok(matrix)
    }
}

fn admit_extent(coatings: usize, elements: usize) -> Result<(), AcousticRealizeError> {
    if coatings > MAX_SHELL_COATINGS
        || coatings
            .checked_mul(elements)
            .is_none_or(|n| n > MAX_COATING_SECTION_ENTRIES)
    {
        return Err(invalid(
            "shell coatings exceed the complete layer/facet coefficient allowance",
        ));
    }
    Ok(())
}

pub(super) fn apply_coatings(
    mut shell: ResolvedShearShell,
    inputs: Vec<ShellCoating>,
    relaxation: Option<relaxation::CoatingRelaxationOptions>,
    gate: &CancelGate,
) -> Result<ResolvedShearShell, AcousticRealizeError> {
    poll(gate)?;
    admit_extent(inputs.len(), shell.mesh.tris.len())?;
    if let Some(options) = relaxation { options.validate()?; }
    shell.coating_relaxation = relaxation;
    if inputs.is_empty() {
        return Ok(shell);
    }
    let mut by_region = vec![Vec::new(); shell.regions.len()];
    let mut records = Vec::with_capacity(inputs.len());
    let mut slices: Vec<CoatingFrequencySlice> = Vec::with_capacity(inputs.len());
    for input in inputs {
        poll(gate)?;
        let r = shell
            .regions
            .iter()
            .position(|r| r.input.region.region.name == input.region)
            .ok_or_else(|| invalid("coating addresses an unknown shell material region"))?;
        if !input.max_storage_relative_drift.is_finite()
            || !(0.0..1.0).contains(&input.max_storage_relative_drift)
        {
            return Err(invalid(
                "coating requires an explicit finite storage-drift bound in [0,1)",
            ));
        }
        let band = input.material.law().band_rad_s();
        let band = (
            band.0.max(shell.regions[r].omega_band.0),
            band.1.min(shell.regions[r].omega_band.1),
        );
        if input.reference_omega_rad_s < band.0 || input.reference_omega_rad_s > band.1 {
            return Err(invalid(
                "coating reference must lie within substrate and film source support",
            ));
        }
        let mut slice = input
            .material
            .at_frequency(input.thickness_m, input.reference_omega_rad_s)
            .map_err(owner)?;
        let (mut storage, loss) = input
            .material
            .law()
            .modulus(input.reference_omega_rad_s)
            .map_err(owner)?;
        if relaxation.is_some() {
            slice.layer = input.material.equilibrium_layer(input.thickness_m).map_err(owner)?;
            storage = input.material.law().equilibrium_modulus_pa();
        }
        let area_m2 = shell.regions[r].area_m2;
        let mass_kg = input.material.density_kg_m3() * input.thickness_m * area_m2;
        if !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(invalid(
                "added coating mass must be positive and representable",
            ));
        }
        by_region[r].push(records.len());
        records.push(ResolvedShellCoating {
            area_m2,
            mass_kg,
            omega_band: band,
            reference_storage_pa: storage,
            reference_viscosity_pa_s: loss / input.reference_omega_rad_s,
            z_bottom_m: 0.0,
            z_top_m: 0.0,
            input,
        });
        slices.push(slice);
    }
    let zero = ShearShellViscosity::new([0.0; 9], [0.0; 9], [0.0; 3]).map_err(owner)?;
    let mut film_viscosity = vec![vec![zero; shell.mesh.tris.len()]; records.len()];
    for (r, indices) in by_region.iter().enumerate() {
        poll(gate)?;
        if indices.is_empty() {
            continue;
        }
        let region = &mut shell.regions[r];
        let h = region.thickness_m;
        let model = region.input.region.model;
        let state = &region.input.region.material;
        let raw = region
            .section
            .principal_matrices()
            .0
            .coefficients()
            .map(|v| v / h);
        let mut gs = if model.is_isotropic() {
            [raw[8]; 2]
        } else {
            [
                plate_property(state, "shear_modulus_13", Pressure::DIMS)?,
                plate_property(state, "shear_modulus_23", Pressure::DIMS)?,
            ]
        };
        if model == PlateMaterialModel::Orthotropic21 {
            gs.swap(0, 1);
        }
        let substrate = LaminaLayer {
            thickness_m: h,
            orientation_rad: 0.0,
            density_kg_m3: region.section.areal_mass_kg_m2() / h,
            material: Some(state.clone()),
            law: LayerLaw::PlaneStress {
                stiffness_pa: raw,
                transverse_shear_pa: gs,
            },
        };
        let mut order: Vec<Option<usize>> = indices
            .iter()
            .rev()
            .copied()
            .filter(|&i| records[i].input.side == CoatingSide::Negative)
            .map(Some)
            .collect();
        let below: f64 = order
            .iter()
            .flatten()
            .map(|&i| records[i].input.thickness_m)
            .sum();
        order.push(None);
        order.extend(
            indices
                .iter()
                .copied()
                .filter(|&i| records[i].input.side == CoatingSide::Positive)
                .map(Some),
        );
        let layers = order
            .iter()
            .map(|i| i.map_or_else(|| substrate.clone(), |i| slices[i].layer.clone()))
            .collect();
        let laminate = LaminateSection::try_new(layers, ReferenceSurface::Offset(below + 0.5 * h))
            .map_err(owner)?;
        // All layers share the substrate frame; films are explicitly isotropic.
        // Preserve the caller's correction, instead of silently substituting
        // the material owner's alternative Whitney correction policy.
        let mut shear = [0.0; 3];
        for geometry in laminate.layer_integrals() {
            let dz = geometry.z_top_m - geometry.z_bottom_m;
            shear[0] += dz * geometry.transverse_shear_pa[0];
            shear[2] += dz * geometry.transverse_shear_pa[1];
        }
        region.section = ShearShellSection::new_coupled(
            SheetMatrix::from_coefficients(*laminate.a_pa_m()).map_err(owner)?,
            *laminate.b_pa_m2(),
            SheetMatrix::from_coefficients(*laminate.d_pa_m3()).map_err(owner)?,
            shear,
            region.input.shear_correction,
            laminate.i0_kg_m2(),
            laminate.i1_kg_m(),
            laminate.i2_kg(),
            0.0,
        )
        .map_err(owner)?;
        region.mass_kg = region.area_m2 * laminate.i0_kg_m2();
        for &i in indices {
            region.omega_band.0 = region.omega_band.0.max(records[i].omega_band.0);
            region.omega_band.1 = region.omega_band.1.min(records[i].omega_band.1);
        }
        if region.omega_band.1 < region.omega_band.0 {
            return Err(invalid("coating stack has no common supported frequency"));
        }
        for (layer_index, &film) in order.iter().enumerate() {
            let Some(i) = film else {
                continue;
            };
            poll(gate)?;
            let geometry = &laminate.layer_integrals()[layer_index];
            if geometry.z_top_m <= geometry.z_bottom_m {
                return Err(invalid(
                    "coating thickness is unresolved at its physical surface offset",
                ));
            }
            records[i].z_bottom_m = geometry.z_bottom_m;
            records[i].z_top_m = geometry.z_top_m;
            let mut values = vec![[0.0; 9]; order.len()];
            values[layer_index] = slices[i].plane_viscosity_pa_s;
            let (a, b, d) = laminate.integrate_layer_matrices(&values).map_err(owner)?;
            let scale = records[i].input.thickness_m * region.input.shear_correction;
            let [gx, gy] = slices[i].transverse_shear_viscosity_pa_s;
            let law = ShearShellViscosity::new_coupled(a, b, d, [scale * gx, 0.0, scale * gy])
                .map_err(owner)?;
            // Isotropic film A/B/D and shear are invariant to in-plane rotation.
            for &e in &region.input.region.region.triangle_indices {
                poll(gate)?;
                film_viscosity[i][e] = law;
            }
        }
        for &e in &region.input.region.region.triangle_indices {
            poll(gate)?;
            let angle = fs_plate::shell::survey::material_axis_angle(
                &shell.mesh.facet(e).map_err(owner)?,
                region.input.axis1_global,
            )
            .map_err(owner)?
                + region.input.region.material_angle_rad;
            shell.sections[e] = region
                .section
                .with_material_angle_rad(angle)
                .map_err(owner)?;
        }
    }
    shell.mass_kg = shell.regions.iter().map(|r| r.mass_kg).sum();
    if !shell.mass_kg.is_finite() || shell.mass_kg <= 0.0 {
        return Err(invalid("coated shell total mass is unrepresentable"));
    }
    shell.coatings = records;
    shell.coating_viscosity = film_viscosity;
    poll(gate)?;
    Ok(shell)
}
