//! Partial-enthalpy transport on the original constituent diffusion transfers.
//!
//! In the barycentric frame the diffusive energy flux is sum(h_k * j_k),
//! with sum(j_k) = 0. A common enthalpy datum therefore carries no heat. The
//! source owner supplies partial SPECIFIC enthalpies in J/kg at each face;
//! molar enthalpy, mixture enthalpy and heat capacity cannot replace them.
//! See Cantera's species/energy equations and their common diffusive fluxes:
//! https://www.cantera.org/3.2/reference/onedim/governing-equations.html.

use super::{
    CancelGate, ChtError, DiffusionObserver, MaterialAdvectionStep, SpeciesDiffusionControls,
    SpeciesDiffusionStep, VoxelDomain, nonnegative, poll, poll_every, refuse,
};

pub mod adaptive;

/// The SAME accepted constituent and thermal inventories after diffusion.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesEnthalpyDiffusionStep<const N: usize> {
    /// Original mass solver's authoritative result and independent mass audits.
    pub species: SpeciesDiffusionStep<N>,
    /// Endpoint total specific enthalpy [J/kg] on the fixed carrier; solids are zero.
    pub specific_enthalpy_j_kg: Vec<f64>,
    /// Signed energy received by each cell through the actual pair transfers [J].
    /// Already included in the endpoint enthalpy; do not deposit it again.
    /// Scalar readback rounds once; internal composition retains its low parts.
    pub diffusive_enthalpy_change_j: Vec<f64>,
    /// Independent physical-unit energy admission.
    pub report: SpeciesEnthalpyDiffusionReport,
}

/// Direct closed-domain thermal storage and transfer admission, in joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeciesEnthalpyDiffusionReport {
    /// Sum of representable m*(h_after-h_before), without subtracting huge totals.
    pub stored_enthalpy_change_j: f64,
    /// Sum of the paired transported-energy increments (zero to roundoff).
    pub transported_enthalpy_change_j: f64,
    /// Rounded nonzero sums of constituent face transfers, weighted by the
    /// largest partial-enthalpy difference. Physical mass-frame diffusion
    /// has exactly zero summed carrier flux.
    pub mass_frame_enthalpy_defect_j: f64,
    /// Absolute local storage/transport mismatches, closed-domain storage and
    /// transport changes, and mass-frame defects. No reference energy scales it.
    pub absolute_enthalpy_defect_j: f64,
}

struct CaloricFace<const N: usize> {
    cells: [usize; 2],
    relative_partial_h: [f64; N],
    relative_partial_h_error: [f64; N],
    maximum_partial_difference: f64,
}

// Neumaier's correction retains a small chemical energy between larger
// opposite signed species contributions, including when the sum changes sign.
#[derive(Clone, Copy, Default)]
struct EnergySum {
    sum: f64,
    correction: f64,
}

impl EnergySum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }

    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

struct EnthalpyTransfers<const N: usize, F> {
    source: F,
    faces: Vec<CaloricFace<N>>,
    increments: Vec<EnergySum>,
    mass_frame_defect: EnergySum,
}

impl<const N: usize, F> DiffusionObserver<N> for EnthalpyTransfers<N, F>
where
    F: Fn([usize; 2], usize) -> Result<[f64; N], ChtError>,
{
    fn storage_bytes(&self, cells: usize, faces: usize) -> Option<usize> {
        // The original owner's mass scratch/output is accounted separately.
        cells
            .checked_mul(core::mem::size_of::<EnergySum>() + 2 * core::mem::size_of::<f64>())?
            .checked_add(faces.checked_mul(core::mem::size_of::<CaloricFace<N>>())?)?
            .checked_add(core::mem::size_of::<Self>()
                + core::mem::size_of::<SpeciesEnthalpyDiffusionStep<N>>())
    }

    fn prepare_storage(&mut self, cells: usize, faces: usize) {
        self.increments.resize(cells, EnergySum::default());
        self.faces.reserve_exact(faces);
    }

    fn prepare(&mut self, cells: [usize; 2], axis: usize) -> Result<(), ChtError> {
        let partial_h = (self.source)(cells, axis)?;
        if partial_h.iter().any(|h| !h.is_finite()) {
            return Err(refuse("partial specific enthalpies must be finite J/kg"));
        }
        // Center before products: a common formation datum cannot amplify
        // floating-point sum(j_k) noise into spurious transported energy.
        let reference = partial_h[0];
        let mut relative_partial_h = [0.0; N];
        let mut relative_partial_h_error = [0.0; N];
        let mut maximum_partial_difference: f64 = 0.0;
        for (k, &h) in partial_h.iter().enumerate() {
            let (difference, tail) = fs_math::eft::two_sum(h, -reference);
            if !difference.is_finite() || !tail.is_finite() {
                return Err(refuse("partial enthalpy differences are unrepresentable"));
            }
            relative_partial_h[k] = difference;
            relative_partial_h_error[k] = tail;
            maximum_partial_difference =
                maximum_partial_difference.max(difference.abs() + tail.abs());
        }
        self.faces.push(CaloricFace {
            cells,
            relative_partial_h,
            relative_partial_h_error,
            maximum_partial_difference,
        });
        Ok(())
    }

    fn transferred(&mut self, face: usize, transfers: [f64; N]) -> Result<(), ChtError> {
        let face = &self.faces[face];
        let mut energy = EnergySum::default();
        let mut carrier_transfer = EnergySum::default();
        for (k, &dm) in transfers.iter().enumerate() {
            // Preserve the original subtraction AND product low parts until
            // different species have cancelled. Summing rounded products alone
            // can report zero heat for a finite change in chemical inventory.
            for dh in [face.relative_partial_h[k], face.relative_partial_h_error[k]] {
                if dh == 0.0 || dm == 0.0 {
                    continue;
                }
                let (work, tail) = fs_math::eft::two_prod(dh, dm);
                if !work.is_finite() || !tail.is_finite() || work == 0.0 {
                    return Err(refuse("partial enthalpy transfer is unrepresentable"));
                }
                energy.add(work);
                energy.add(tail);
            }
            carrier_transfer.add(dm);
        }
        if !energy.value().is_finite() {
            return Err(refuse("diffusive enthalpy sum overflowed"));
        }
        let frame_defect = carrier_transfer.value().abs() * face.maximum_partial_difference;
        if !frame_defect.is_finite() {
            return Err(refuse("diffusive mass-frame energy defect overflowed"));
        }
        self.mass_frame_defect.add(frame_defect);
        for (cell, sign) in [(face.cells[0], -1.0), (face.cells[1], 1.0)] {
            let increment = &mut self.increments[cell];
            // A later face or sweep may cancel the large part. Keep its
            // remainder until all transfers into this cell have met.
            increment.add(sign * energy.sum);
            increment.add(sign * energy.correction);
            if !increment.value().is_finite() {
                return Err(refuse("diffusive enthalpy inventory overflowed"));
            }
        }
        Ok(())
    }
}

impl MaterialAdvectionStep {
    /// Diffuse the original species AND their partial enthalpy on one candidate.
    ///
    /// The initial specific enthalpy is total caloric storage on this accepted
    /// carrier, with the SAME source/formation convention as the partial
    /// enthalpies. The callback runs once per active positive-axis fluid-fluid
    /// face with its ordered cell indices and axis (0, 1, or 2). It supplies
    /// source-evaluated face partial specific enthalpies [J/kg], frozen for
    /// this interval. No composition, Cp, temperature or material identity is
    /// inferred from species indices. A common datum cancels before products.
    /// Formation and latent contributions already in h_k must not be credited
    /// again as heat sources, nor may pressure-flow work be added twice.
    ///
    /// The original positive mass relaxation, half-cell conductances, ordered
    /// sweeps and constituent audits are reused without another mass solve.
    /// Each actual signed species transfer supplies sum(h_k*dm_k), debited and
    /// credited once on the two original cells. Subtraction/product low parts
    /// survive species cancellation. Compensated increments precede
    /// a centered m*h update. Unrepresentable storage work and rounded
    /// mass-frame defects are charged to the absolute joule budget, never
    /// hidden by a relative formation-energy scale. The original storage
    /// budget includes thermal scratch, face coefficients and outputs.
    ///
    /// This is a frozen partial-enthalpy diffusion component, not an EOS or
    /// a coupled nonisothermal integrator. The caller must invert the resulting
    /// total enthalpy WITH the resulting composition using its original law,
    /// refine/re-evaluate face properties as that state changes, and commit
    /// both inventories together. No temperature, reaction, pressure-work,
    /// Soret, separate heat-of-mixing or second-order coupled-time claim is made.
    /// Failure/cancellation publishes neither inventory and mutates no input.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_species_with_enthalpy<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        face_partial_enthalpy: impl Fn([usize; 2], usize) -> Result<[f64; N], ChtError>,
        controls: SpeciesDiffusionControls<N>,
        maximum_enthalpy_defect_j: f64,
        gate: &CancelGate,
    ) -> Result<SpeciesEnthalpyDiffusionStep<N>, ChtError> {
        self.diffuse_species_with_enthalpy_accounted(
            domain, species_mass_kg, specific_enthalpy_j_kg, diffusivity_m2_s,
            duration_s, face_partial_enthalpy, controls, maximum_enthalpy_defect_j, gate,
        ).map(|(step, _)| step)
    }

    // Internal composition retains the original transported-energy low parts.
    // The public scalar readback rounds only at its output boundary; an outer
    // physical integrator must not reconstruct these increments from it.
    #[allow(clippy::too_many_arguments)]
    fn diffuse_species_with_enthalpy_accounted<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        diffusivity_m2_s: &[f64],
        duration_s: f64,
        face_partial_enthalpy: impl Fn([usize; 2], usize) -> Result<[f64; N], ChtError>,
        controls: SpeciesDiffusionControls<N>,
        maximum_enthalpy_defect_j: f64,
        gate: &CancelGate,
    ) -> Result<(SpeciesEnthalpyDiffusionStep<N>, Vec<EnergySum>), ChtError> {
        poll(gate)?;
        nonnegative(maximum_enthalpy_defect_j)?;
        if specific_enthalpy_j_kg.len() != domain.cell_count() {
            return Err(refuse("enthalpy field must match the diffusion domain"));
        }
        for (cell, &h) in specific_enthalpy_j_kg.iter().enumerate() {
            poll_every(gate, cell)?;
            if !h.is_finite() || (!domain.is_fluid(cell) && h != 0.0) {
                return Err(refuse("fluid enthalpy must be finite and solid entries zero"));
            }
        }
        let mut thermal = EnthalpyTransfers {
            source: face_partial_enthalpy,
            faces: Vec::new(),
            increments: Vec::new(),
            mass_frame_defect: EnergySum::default(),
        };
        let species = self.diffuse_species_observed(
            domain, species_mass_kg, diffusivity_m2_s, duration_s, controls, gate,
            |_| poll(gate), &mut thermal,
        )?;
        let mut specific_h = Vec::with_capacity(specific_enthalpy_j_kg.len());
        let mut changes = Vec::with_capacity(specific_enthalpy_j_kg.len());
        let mut stored = EnergySum::default();
        let mut transported = EnergySum::default();
        let mut defect = EnergySum::default();
        for (cell, &initial_h) in specific_enthalpy_j_kg.iter().enumerate() {
            poll_every(gate, cell)?;
            let increment = thermal.increments[cell];
            let exchange = increment.value();
            let mass = self.mass_kg[cell];
            let h = if domain.is_fluid(cell) {
                let mut candidate = EnergySum::default();
                candidate.add(initial_h);
                candidate.add(increment.sum / mass);
                candidate.add(increment.correction / mass);
                candidate.value()
            } else {
                initial_h
            };
            // Audit the represented endpoint with retained subtraction and
            // product remainders too. Rounded m*delta_h can equal the requested
            // heat while the actual represented inventory misses whole joules.
            let (delta_h, delta_h_error) = fs_math::eft::two_sum(h, -initial_h);
            let (actual, actual_error) = fs_math::eft::two_prod(mass, delta_h);
            let (extra, extra_error) = fs_math::eft::two_prod(mass, delta_h_error);
            let represented = [actual, actual_error, extra, extra_error];
            if !h.is_finite() || represented.iter().any(|term| !term.is_finite()) {
                return Err(refuse("diffusive enthalpy storage is unrepresentable"));
            }
            let mut local_defect = EnergySum::default();
            for term in represented {
                stored.add(term);
                local_defect.add(term);
            }
            local_defect.add(-increment.sum);
            local_defect.add(-increment.correction);
            specific_h.push(h);
            changes.push(exchange);
            transported.add(increment.sum);
            transported.add(increment.correction);
            defect.add(local_defect.value().abs());
        }
        let absolute_enthalpy_defect_j = defect.value() + stored.value().abs()
            + transported.value().abs() + thermal.mass_frame_defect.value();
        nonnegative(absolute_enthalpy_defect_j)?;
        if absolute_enthalpy_defect_j > maximum_enthalpy_defect_j {
            return Err(refuse("diffusion exceeds its enthalpy balance budget"));
        }
        poll(gate)?;
        let step = SpeciesEnthalpyDiffusionStep {
            species,
            specific_enthalpy_j_kg: specific_h,
            diffusive_enthalpy_change_j: changes,
            report: SpeciesEnthalpyDiffusionReport {
                stored_enthalpy_change_j: stored.value(),
                transported_enthalpy_change_j: transported.value(),
                mass_frame_enthalpy_defect_j: thermal.mass_frame_defect.value(),
                absolute_enthalpy_defect_j,
            },
        };
        Ok((step, thermal.increments))
    }
}
