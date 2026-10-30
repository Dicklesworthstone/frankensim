//! Gas heat diffusion coupled to the original reciprocal compression field.
//!
//! For expansion coefficients e_j=xi_j/sqrt(Lambda_j), let H_lj be the
//! L2 overlap of the normalized pressure shape with Dirichlet heat mode l.
//! With B=gamma*p0 and D=B*(1-1/gamma), the exact finite-basis constitutive law is
//!   p_hat = -B e + sqrt(D) H^T z,
//!   z'_l = a_l (sqrt(D) H_l e-z_l), a_l=kappa*lambda_l/(rho*cv).
//! Its storage is e^T(B I-D H^T H)e/2 + sum_l(sqrt(D)H_l e-z_l)^2/2.
//! Thus all missing heat modes remain adiabatic, never silently isothermal.
//! Their compression-projector remainder is explicitly bounded before binding.
//! Existing Maxwell states, Cholesky and the original network own all solves.
//!
//! Reference equations: linearized ideal-gas energy and EOS in COMSOL's
//! Acoustic Perturbation and Linearization documentation, equations 6-6..6-10:
//! https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_thermo.09.79.html
//! This is a direct Galerkin reduction of those equations, not a frequency-fit
//! or a modal-Q assignment. Solid walls are fixed-temperature reservoirs and
//! openings have zero diffusive heat flux. No viscous wall shear, solid thermal
//! impedance, mean flow, aperture enthalpy transport or finite heating follows.

use super::*;
use crate::render::plate::impact::cavity::network::ModalCavityParts;
use fs_feec::acoustics::helmholtz::diffusion::DiffusionModes;
use fs_phs::RelaxationBranch;

mod molecular;
pub use molecular::ShellCavityMolecularOptions;

/// Thermal FE/window admission, independent of instrument identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellCavityThermalOptions {
    /// Keep the entire positive heat spectrum through this decay rate [1/s].
    pub maximum_decay_rate_per_s: f64,
    /// Maximum complete heat eigenwindow and added memory count, 1..=256.
    pub maximum_thermal_modes: usize,
    /// Maximum original volume-space DOFs, 1..=8192.
    pub maximum_thermal_dofs: usize,
    /// Sparse products and dense projection/factor work, at most 10^9 terms.
    pub maximum_projection_terms: usize,
    /// Bound on the unresolved normalized compression projector, strictly <1.
    /// This measures FE+modal representation, not physical instrument error.
    pub maximum_unresolved_fraction: f64,
}
impl ShellCavityThermalOptions {
    pub fn validate(&self) -> Result<(), ShellCavityError> {
        if !self.maximum_decay_rate_per_s.is_finite() || self.maximum_decay_rate_per_s <= 0.0
            || !(1..=256).contains(&self.maximum_thermal_modes)
            || !(1..=8192).contains(&self.maximum_thermal_dofs)
            || !(1..=1_000_000_000).contains(&self.maximum_projection_terms)
            || !self.maximum_unresolved_fraction.is_finite()
            || self.maximum_unresolved_fraction <= 0.0 || self.maximum_unresolved_fraction >= 1.0
        { return Err(bad("thermal cavity requires positive bounded spectrum, memory, work and remainder declarations")); }
        Ok(())
    }
}

/// Retained original gas diffusion and work-conjugate pressure projections.
#[derive(Clone, Debug)]
pub struct ShellCavityThermalLoss {
    options: ShellCavityThermalOptions,
    diffusion: DiffusionModes,
    diffusivity_m2_s: f64,
    unresolved_fraction: f64,
    /// One physical decay rate for each active compression memory.
    rates: Vec<f64>,
    mode_indices: Vec<usize>,
    /// Rows multiplying original mass-normalized mechanical q, in 1/s.
    memory_rows: Vec<Vec<f64>>,
    /// Original normalized pressure projector, in full diffusion mode order.
    overlaps: Vec<Vec<f64>>,
    /// Actual pressure [Pa] per energy-normalized thermal memory.
    pressure_rows: Vec<Vec<f64>>,
    /// Cholesky rows of the remaining equilibrium compression energy.
    equilibrium_rows: Vec<Vec<f64>>,
    /// Original compression-volume normalizations for physical constitutive observation.
    pressure_norms: Vec<f64>,
    instantaneous_bulk_modulus_pa: f64,
    instantaneous_pressure_scale: f64,
    molecular: Option<fs_material::relaxation::RelaxationState>,
}
impl ShellCavityThermalLoss {
    pub const fn options(&self) -> &ShellCavityThermalOptions { &self.options }
    pub const fn diffusion_modes(&self) -> &DiffusionModes { &self.diffusion }
    pub const fn diffusivity_m2_s(&self) -> f64 { self.diffusivity_m2_s }
    pub const fn unresolved_fraction_bound(&self) -> f64 { self.unresolved_fraction }
    pub fn decay_rates_per_s(&self) -> &[f64] { &self.rates }
    /// Original eigenmode for each wall-diffusion memory, in memory-prefix order.
    /// Joint caloric branches can repeat an index. Molecular-only histories in
    /// the unresolved spatial complement follow this prefix and have no heat
    /// eigenmode. Exactly zero projections remain in diffusion_modes().
    pub fn memory_mode_indices(&self) -> &[usize] { &self.mode_indices }
    pub fn memory_states(&self) -> usize { self.memory_rows.len() }
    pub fn normalized_pressure_overlaps(&self) -> &[Vec<f64>] { &self.overlaps }
    /// Exact source gas and internal reservoirs used by the joint caloric law.
    pub const fn molecular_relaxation(&self) -> Option<&fs_material::relaxation::RelaxationState> {
        self.molecular.as_ref()
    }
    /// Molecular histories in spatial directions not resolved by the heat basis.
    pub fn molecular_complement_memory_states(&self) -> usize {
        self.memory_rows.len() - self.mode_indices.len()
    }
    /// Instantaneous compression with wall diffusion and internal reservoirs frozen [Pa].
    pub const fn instantaneous_bulk_modulus_pa(&self) -> f64 {
        self.instantaneous_bulk_modulus_pa
    }
    /// Constitutive bulk matrix in the original mass-normalized pressure basis,
    /// using exp(-i omega t). This observes the SAME storage and relaxation rows
    /// as the native mechanical runtime; it is not a separate cavity field solve.
    /// Pressure-basis and geometric convergence remain the original FE obligations.
    pub fn normalized_bulk_response(&self, omega_rad_s: f64)
        -> Result<Vec<fs_math::c64::C64>, ShellCavityError>
    {
        use fs_math::c64::C64;
        if !omega_rad_s.is_finite() || omega_rad_s < 0.0 {
            return Err(bad("caloric response requires finite nonnegative frequency"));
        }
        let n = self.pressure_norms.len();
        let mut result = vec![C64::ZERO; n*n];
        for i in 0..n { result[i*n+i] = C64::from_re(self.instantaneous_bulk_modulus_pa); }
        for (row, &rate) in self.pressure_rows.iter().zip(&self.rates) {
            let response = C64::from_re(rate) / C64::new(rate, -omega_rad_s);
            for i in 0..n { for j in 0..n {
                result[i*n+j] = result[i*n+j] - response.scale(
                    row[i]*self.pressure_norms[i]*row[j]*self.pressure_norms[j]);
            }}
        }
        if result.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
            return Err(bad("caloric response is unrepresentable"));
        }
        Ok(result)
    }

    pub(super) fn prepare(
        fe: &HelmholtzFe<'_>, gas: &GasState, shapes: &[Vec<f64>],
        coupling: &CavityCoupling, openings: &[openings::PreparedOpening],
        options: ShellCavityThermalOptions, eigensolver: &SliceOptions, gate: &CancelGate,
    ) -> Result<Self, ShellCavityError> {
        check_cancel(gate)?; options.validate()?;
        gas.validate_coherence().map_err(|_|bad("thermal cavity requires its original coherent gas"))?;
        let pressure_count = coupling.cavity_modes();
        if fe.ndof() > options.maximum_thermal_dofs || shapes.len() != pressure_count {
            return Err(bad("thermal cavity lost the original pressure basis or exceeded the FE allowance"));
        }
        // Include every sparse product, pressure-pair matrix and actual
        // original mechanical row before calling the common eigensolver.
        let nonzeros = (0..fe.ndof()).try_fold(0usize, |n,row| n.checked_add(fe.mass().row(row).0.len()));
        let heat = options.maximum_thermal_modes;
        let work = nonzeros.and_then(|n| n.checked_mul(heat.checked_mul(heat.checked_add(pressure_count)?)?))
            .and_then(|n| pressure_count.checked_mul(pressure_count)
                .and_then(|p| p.checked_mul(heat.checked_add(pressure_count)?.checked_add(coupling.total_modes())?))
                .and_then(|p| n.checked_add(p)));
        if work.is_none_or(|n| n > options.maximum_projection_terms) {
            return Err(bad("complete heat-mode projection exceeds its original work allowance"));
        }
        let cv = gas.specific_heat_cp / gas.gamma;
        let diffusivity = gas.thermal_conductivity / (gas.density * cv);
        let cutoff = options.maximum_decay_rate_per_s / diffusivity;
        if !diffusivity.is_finite() || diffusivity <= 0.0 || !cutoff.is_finite() || cutoff <= 0.0 {
            return Err(bad("gas heat diffusivity or complete eigenwindow is unrepresentable"));
        }
        let mut aperture_faces = std::collections::BTreeSet::new();
        for opening in openings { for &face in &opening.source.boundary_faces {
            let mut canonical = face; canonical.sort_unstable(); aperture_faces.insert(canonical);
        }}
        let walls: Vec<_> = fe.boundary_faces().into_iter()
            .filter(|face| !aperture_faces.contains(face)).collect();
        let diffusion = fe.isothermal_diffusion_modes(&walls, cutoff, heat, eigensolver,
            ||gate.is_requested()).map_err(ShellCavityError::Acoustic)?;
        check_cancel(gate)?;
        let storage = coupling.compression_storage();
        let norms: Vec<_> = storage.iter().map(|spring|det::sqrt(spring.volume_m3)).collect();
        let mut overlaps = Vec::with_capacity(diffusion.shapes.len());
        for (index, thermal) in diffusion.shapes.iter().enumerate() {
            for (other, shape) in diffusion.shapes[..=index].iter().enumerate() {
                let gram = mass_product(fe,thermal,shape,gate)?;
                if (gram-if other==index {1.0} else {0.0}).abs() > 1e-7 {
                    return Err(bad("thermal eigenvectors must retain the original mass-orthonormal basis"));
                }
            }
            let mut row = Vec::with_capacity(pressure_count);
            for (shape, norm) in shapes.iter().zip(&norms) {
                row.push(mass_product(fe,thermal,shape,gate)? / norm);
            }
            if row.iter().any(|v|!v.is_finite()) { return Err(bad("thermal compression projection overflow")); }
            overlaps.push(row);
        }
        let bulk = gas.density * gas.sound_speed * gas.sound_speed;
        let fraction = 1.0 - 1.0 / gas.gamma;
        let increment = bulk * fraction;
        let root_bulk = det::sqrt(bulk);
        let root_increment = det::sqrt(increment);
        let mut equilibrium = vec![0.0; pressure_count*pressure_count];
        let mut unresolved = 0.0_f64;
        for i in 0..pressure_count {
            check_cancel(gate)?;
            let mut row_sum = 0.0;
            for j in 0..pressure_count {
                let projected: f64 = overlaps.iter().map(|row|row[i]*row[j]).sum();
                let identity = if i==j {1.0} else {0.0};
                row_sum += (identity-projected).abs();
                equilibrium[i*pressure_count+j] = identity - fraction*projected;
            }
            unresolved = unresolved.max(row_sum);
        }
        if !unresolved.is_finite() || unresolved > options.maximum_unresolved_fraction {
            return Err(bad("thermal FE and eigenwindow leave excessive unresolved compression; refine the actual volume or retain more heat modes"));
        }
        let factor = fs_la::factor::cholesky(&equilibrium,pressure_count)
            .map_err(|_|bad("original projected isothermal compression form is not positive definite"))?;
        let mut equilibrium_rows = vec![vec![0.0;coupling.total_modes()];pressure_count];
        for (row, target) in equilibrium_rows.iter_mut().enumerate() {
            check_cancel(gate)?;
            for (j,spring) in storage.iter().enumerate() {
                let coefficient = root_bulk * factor.l(j,row) / norms[j];
                for (target,&shape) in target.iter_mut().zip(&spring.areas) { *target += coefficient*shape; }
            }
        }
        let mut rates = Vec::new(); let mut mode_indices = Vec::new();
        let mut memory_rows = Vec::new(); let mut pressure_rows = Vec::new();
        for (index, overlap) in overlaps.iter().enumerate() {
            check_cancel(gate)?;
            let mut row = vec![0.0;coupling.total_modes()];
            let mut pressure = Vec::with_capacity(pressure_count);
            for (j,spring) in storage.iter().enumerate() {
                let coefficient = root_increment * overlap[j] / norms[j];
                pressure.push(coefficient);
                for (value,&shape) in row.iter_mut().zip(&spring.areas) { *value += coefficient*shape; }
            }
            if row.iter().chain(&pressure).any(|v|!v.is_finite()) {
                return Err(bad("thermal memory or physical pressure projection overflow"));
            }
            if row.iter().all(|v|*v==0.0) { continue; }
            let rate = diffusivity*diffusion.eigenvalues_per_m2[index];
            if !rate.is_finite() || rate <= 0.0 || !(1.0/rate).is_finite() || 1.0/rate <= 0.0 {
                return Err(bad("thermal decay time is unrepresentable"));
            }
            rates.push(rate); mode_indices.push(index); memory_rows.push(row); pressure_rows.push(pressure);
        }
        if memory_rows.is_empty() || equilibrium_rows.iter().flatten().any(|v|!v.is_finite()) {
            return Err(bad("thermal source must retain representable compression memory and equilibrium storage"));
        }
        check_cancel(gate)?;
        Ok(Self { options, diffusion, diffusivity_m2_s:diffusivity, unresolved_fraction:unresolved,
            rates, mode_indices, memory_rows, overlaps, pressure_rows, equilibrium_rows,
            pressure_norms:norms, instantaneous_bulk_modulus_pa:bulk,
            instantaneous_pressure_scale:1.0, molecular:None })
    }

    pub(super) fn replace_compression_links(&self, image: &mut ModalCavityParts,
        coupling: &CavityCoupling, gate: &CancelGate) -> Result<(),ShellCavityError>
    {
        let count = coupling.cavity_modes(); let structural = coupling.structural_modes();
        let extra = coupling.total_modes()-structural;
        if image.connections.len() < count || self.equilibrium_rows.len()!=count {
            return Err(bad("thermal lowering lost original compression links"));
        }
        for (link,row) in image.connections[..count].iter_mut().zip(&self.equilibrium_rows) {
            check_cancel(gate)?;
            let right = if extra==0 {
                ModalAttachment { component:0, shapes:vec![0.0;structural] }
            } else {
                ModalAttachment { component:1, shapes:row[structural..].iter().map(|v|-v).collect() }
            };
            *link = ModalConnection { left:ModalAttachment { component:0, shapes:row[..structural].to_vec() },
                right, stiffness_n_m:1.0, damping_n_s_m:0.0, rest_extension_m:0.0 };
        }
        Ok(())
    }

    pub(super) fn correct_pressure(&self, memory: &[f64], pressure: &mut [f64])
        -> Result<(),ShellCavityError>
    {
        let start = memory.len().checked_sub(self.memory_states())
            .ok_or_else(||bad("cavity pressure lost its accepted thermal history"))?;
        if pressure.len()!=self.equilibrium_rows.len() || memory[start..].iter().any(|z|!z.is_finite()) {
            return Err(bad("cavity pressure requires complete finite thermal and pressure coordinates"));
        }
        for p in pressure.iter_mut() { *p *= self.instantaneous_pressure_scale; }
        for (row,&z) in self.pressure_rows.iter().zip(&memory[start..]) {
            for (p,coefficient) in pressure.iter_mut().zip(row) { *p += coefficient*z; }
        }
        if pressure.iter().any(|v|!v.is_finite()) { return Err(bad("thermal cavity pressure overflow")); }
        Ok(())
    }
}

impl MaterialShellCavity {
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn thermal_memory_branches(
        &self, source: &MaterialStringBankSource, system: &FrictionModalSystem, gate: &CancelGate,
    ) -> Result<Vec<RelaxationBranch>,ReciprocalStringError> {
        let Some(thermal) = &self.thermal_loss else { return Ok(Vec::new()); };
        let body = source.strings.len();
        let count = 1+usize::from(self.coupling.total_modes()>self.coupling.structural_modes());
        let components = system.components().get(body..body+count)
            .ok_or_else(||invalid("thermal history lost its original body and cavity components"))?;
        if components[0].modes().len()!=self.coupling.structural_modes()
            || components.iter().map(|c|c.modes().len()).sum::<usize>()!=self.coupling.total_modes() {
            return Err(invalid("thermal history cannot be remapped onto a different cavity basis"));
        }
        let start:usize = system.components()[..body].iter().map(|c|c.modes().len()).sum();
        let dimension = system.mode_count().checked_mul(2)
            .ok_or_else(||invalid("complete thermal history dimension overflow"))?;
        if thermal.memory_states().checked_mul(dimension)
            .is_none_or(|n|n>source.options.network.max_setup_terms) {
            return Err(invalid("complete thermal history exceeds the original network setup allowance"));
        }
        let mut branches = Vec::with_capacity(thermal.memory_states());
        for (row,&rate) in thermal.memory_rows.iter().zip(&thermal.rates) {
            poll(gate)?;
            let mut projection = vec![0.0;dimension];
            for (index,&value) in row.iter().enumerate() { projection[2*(start+index)]=value; }
            branches.push(RelaxationBranch { projection,stiffness:1.0,relaxation_time_s:1.0/rate });
        }
        Ok(branches)
    }
}

#[cfg(test)]
mod tests;
