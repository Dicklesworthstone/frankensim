//! Partial enthalpy on the original fixed-composition reservoir maps.
//!
//! Every actual boundary constituent transfer is observed once, positive INTO
//! the fluid cell. The face law supplies partial specific enthalpies in the
//! same energy convention as cell storage. This is diffusive exchange on a
//! fixed carrier, not bulk inflow or a second advective enthalpy deposit.

pub mod adaptive;

use super::super::boundary::{
    ReservoirObserver, ReservoirTransfers, SpeciesBoundaryDiffusionReport, SpeciesDiffusionBoundary,
};
use super::{
    CancelGate, ChtError, DiffusionObserver, EnergySum, EnthalpyTransfers, MaterialAdvectionStep,
    PartialEnthalpy, SpeciesDiffusionControls, SpeciesEnthalpyDiffusionStep, VoxelDomain,
    nonnegative, poll, poll_every, refuse, validate_enthalpy,
};

/// Energy received through one declared reservoir face, in input order.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesBoundaryEnthalpyTransfer {
    pub cell: usize,
    pub face: usize,
    /// Signed actual diffusive heat received by the cell [J]. This is already
    /// present in its endpoint enthalpy; do not add it a second time.
    pub enthalpy_into_cell_j: f64,
    /// Sum of positive NET map heats [J], after constituent cancellation.
    /// This is not a sum of reference-dependent individual species energies.
    pub positive_enthalpy_into_cell_j: f64,
    /// Magnitude of negative NET map heats [J]. Both gross signs may be nonzero
    /// when an interior exchange or another reservoir reverses the face flux.
    pub negative_enthalpy_into_cell_j: f64,
}

/// Independent external energy account; scalar readbacks round only here.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesBoundaryEnthalpyDiffusionReport {
    /// Includes inactive declared faces with exactly zero exchanged heat.
    pub faces: Vec<SpeciesBoundaryEnthalpyTransfer>,
    pub enthalpy_into_domain_j: f64,
    pub positive_enthalpy_into_domain_j: f64,
    pub negative_enthalpy_into_domain_j: f64,
}

/// One atomic species/enthalpy candidate with its original reservoir receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct BoundarySpeciesEnthalpyDiffusionStep<const N: usize> {
    pub step: SpeciesEnthalpyDiffusionStep<N>,
    pub boundary: SpeciesBoundaryDiffusionReport<N>,
    pub thermal_boundary: SpeciesBoundaryEnthalpyDiffusionReport,
}

// Internal physical composition must retain all original transfer low parts.
// In particular, neither a large cancelling formation energy nor a large
// outward species transfer may erase a smaller later inward transfer.
#[derive(Clone, Copy)]
pub(super) struct BoundaryEnthalpySums<const N: usize> {
    pub(super) signed: EnergySum,
    pub(super) positive: EnergySum,
    pub(super) negative: EnergySum,
    pub(super) mass_signed: [EnergySum; N],
    pub(super) mass_incoming: [EnergySum; N],
    pub(super) mass_outgoing: [EnergySum; N],
}

impl<const N: usize> Default for BoundaryEnthalpySums<N> {
    fn default() -> Self {
        Self {
            signed: EnergySum::default(),
            positive: EnergySum::default(),
            negative: EnergySum::default(),
            mass_signed: [EnergySum::default(); N],
            mass_incoming: [EnergySum::default(); N],
            mass_outgoing: [EnergySum::default(); N],
        }
    }
}

pub(super) struct BoundaryEnthalpyAccount<const N: usize> {
    pub(super) increments: Vec<EnergySum>,
    pub(super) faces: Vec<BoundaryEnthalpySums<N>>,
    pub(super) external: EnergySum,
    pub(super) positive: EnergySum,
    pub(super) negative: EnergySum,
    pub(super) mass_into_domain: [EnergySum; N],
}

/// Inspect the already-admitted canonical boundary list without allocating
/// topology or re-querying constitutive laws. Adjacent entries suffice for
/// shared-cell reservoirs because the original owner requires (cell, face)
/// order. A positive-D fluid neighbor supplies a competing interior map.
pub(super) fn boundary_maps_compete<const N: usize>(
    domain: &VoxelDomain,
    diffusivity: &[f64],
    boundaries: &[SpeciesDiffusionBoundary<N>],
    gate: &CancelGate,
) -> Result<bool, ChtError> {
    let mut previous_active = None;
    for (index, boundary) in boundaries.iter().enumerate() {
        poll_every(gate, index)?;
        let cell = boundary.cell;
        if diffusivity[cell] == 0.0 {
            continue;
        }
        if previous_active == Some(cell) {
            return Ok(true);
        }
        previous_active = Some(cell);
        for face in 0..6 {
            if domain
                .neighbor(cell, face)
                .is_some_and(|n| domain.is_fluid(n) && diffusivity[n] > 0.0)
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

struct BoundaryCaloricFace<const N: usize> {
    input_index: usize,
    cell: usize,
    partials: PartialEnthalpy<N>,
}

struct BoundaryEnthalpyTransfers<const N: usize, F, B> {
    thermal: EnthalpyTransfers<N, F>,
    source: B,
    faces: Vec<BoundaryCaloricFace<N>>,
    account: BoundaryEnthalpyAccount<N>,
}

impl<const N: usize, F, B> DiffusionObserver<N> for BoundaryEnthalpyTransfers<N, F, B>
where
    F: Fn([usize; 2], usize) -> Result<[f64; N], ChtError>,
    B: Fn(usize) -> Result<[f64; N], ChtError>,
{
    fn storage_bytes(&self, cells: usize, faces: usize) -> Option<usize> {
        self.thermal.storage_bytes(cells, faces)
    }
    fn prepare_storage(&mut self, cells: usize, faces: usize) {
        self.thermal.prepare_storage(cells, faces);
    }
    fn prepare(&mut self, cells: [usize; 2], axis: usize) -> Result<(), ChtError> {
        self.thermal.prepare(cells, axis)
    }
    fn transferred(&mut self, face: usize, transfers: [f64; N]) -> Result<(), ChtError> {
        self.thermal.transferred(face, transfers)
    }
}

impl<const N: usize, F, B> ReservoirObserver<N> for BoundaryEnthalpyTransfers<N, F, B>
where
    F: Fn([usize; 2], usize) -> Result<[f64; N], ChtError>,
    B: Fn(usize) -> Result<[f64; N], ChtError>,
{
    fn reservoir_storage_bytes(&self, declared: usize, active: usize) -> Option<usize> {
        active
            .checked_mul(core::mem::size_of::<BoundaryCaloricFace<N>>())?
            .checked_add(declared.checked_mul(
                core::mem::size_of::<BoundaryEnthalpySums<N>>()
                    + core::mem::size_of::<SpeciesBoundaryEnthalpyTransfer>(),
            )?)?
            .checked_add(core::mem::size_of::<Self>())?
            .checked_add(core::mem::size_of::<BoundarySpeciesEnthalpyDiffusionStep<N>>())
    }

    fn prepare_reservoir_storage(&mut self, declared: usize, active: usize) {
        self.faces.reserve_exact(active);
        self.account
            .faces
            .resize_with(declared, BoundaryEnthalpySums::default);
    }

    fn prepare_reservoir(
        &mut self,
        input_index: usize,
        boundary: &SpeciesDiffusionBoundary<N>,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        let partials = (self.source)(input_index)?;
        poll(gate)?;
        self.faces.push(BoundaryCaloricFace {
            input_index,
            cell: boundary.cell,
            partials: PartialEnthalpy::new(partials)?,
        });
        Ok(())
    }

    fn reservoir_transferred(
        &mut self,
        active_index: usize,
        transfers: [f64; N],
    ) -> Result<(), ChtError> {
        let face = &self.faces[active_index];
        let (energy, frame_defect) = face.partials.transfer(transfers)?;
        self.thermal.mass_frame_defect.add(frame_defect);
        add_retained(&mut self.thermal.increments[face.cell], energy)?;
        add_retained(&mut self.account.external, energy)?;
        let receipt = &mut self.account.faces[face.input_index];
        add_retained(&mut receipt.signed, energy)?;
        if energy.value() >= 0.0 {
            add_retained(&mut receipt.positive, energy)?;
            add_retained(&mut self.account.positive, energy)?;
        } else {
            let magnitude = EnergySum {
                sum: -energy.sum,
                correction: -energy.correction,
            };
            add_retained(&mut receipt.negative, magnitude)?;
            add_retained(&mut self.account.negative, magnitude)?;
        }
        for (k, &dm) in transfers.iter().enumerate() {
            receipt.mass_signed[k].add(dm);
            self.account.mass_into_domain[k].add(dm);
            if dm >= 0.0 {
                receipt.mass_incoming[k].add(dm);
            } else {
                receipt.mass_outgoing[k].add(-dm);
            }
            for sum in [
                receipt.mass_signed[k],
                receipt.mass_incoming[k],
                receipt.mass_outgoing[k],
                self.account.mass_into_domain[k],
            ] {
                if !sum.value().is_finite() {
                    return Err(refuse("retained reservoir mass transfer overflowed"));
                }
            }
        }
        Ok(())
    }
}

impl<const N: usize, F, B> BoundaryEnthalpyTransfers<N, F, B> {
    fn report(
        &self,
        boundaries: &[SpeciesDiffusionBoundary<N>],
        gate: &CancelGate,
    ) -> Result<SpeciesBoundaryEnthalpyDiffusionReport, ChtError> {
        let mut faces = Vec::with_capacity(boundaries.len());
        for (index, (boundary, sums)) in boundaries.iter().zip(&self.account.faces).enumerate() {
            poll_every(gate, index)?;
            nonnegative(sums.positive.value())?;
            nonnegative(sums.negative.value())?;
            faces.push(SpeciesBoundaryEnthalpyTransfer {
                cell: boundary.cell,
                face: boundary.face,
                enthalpy_into_cell_j: sums.signed.value(),
                positive_enthalpy_into_cell_j: sums.positive.value(),
                negative_enthalpy_into_cell_j: sums.negative.value(),
            });
        }
        poll(gate)?;
        Ok(SpeciesBoundaryEnthalpyDiffusionReport {
            faces,
            enthalpy_into_domain_j: self.account.external.value(),
            positive_enthalpy_into_domain_j: self.account.positive.value(),
            negative_enthalpy_into_domain_j: self.account.negative.value(),
        })
    }
}

fn add_retained(target: &mut EnergySum, value: EnergySum) -> Result<(), ChtError> {
    target.add(value.sum);
    target.add(value.correction);
    if !target.value().is_finite() {
        return Err(refuse("reservoir enthalpy transfer overflowed"));
    }
    Ok(())
}

impl MaterialAdvectionStep {
    /// Diffuse species and total enthalpy through original interior/reservoir maps.
    ///
    /// The interior callback supplies partial specific enthalpies [J/kg] once
    /// per active positive-axis fluid face, as in `diffuse_species_with_enthalpy`.
    /// The boundary callback receives the canonical index into `boundaries`,
    /// once per active declared reservoir. It supplies ONE physical face law
    /// shared by all counterdiffusing species, in the same constituent order
    /// and energy convention as cell storage. Zero-D faces are never queried.
    /// Film resistance and half-cell conductance remain with the original mass
    /// owner; this entry infers neither a face temperature nor a heat-film law.
    ///
    /// Each actual boundary dm is positive into the cell. Centered partial
    /// differences, retained subtraction/product remainders and compensated
    /// increments preserve heat across large chemical-energy cancellation.
    /// The original endpoint audit now compares both storage and transported
    /// energy against independently accumulated EXTERNAL heat, retaining its
    /// low parts; local represented-storage and barycentric defects still cost
    /// absolute joules. Boundary heat is already included in the returned h.
    /// Gross signs describe net map heat after species cancellation, so they
    /// remain invariant under a common partial-enthalpy datum. Reported scalar
    /// heats round once; internal physical composition retains their low parts.
    ///
    /// All original mass, resource, geometry and cancellation gates apply.
    /// Storage admission includes thermal face caches and both retained and
    /// public receipts before allocating them. No inputs are mutated; a failure
    /// publishes neither inventory. Sources are frozen for this interval: the
    /// caller owns temperature/composition inversion and property refresh.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_species_with_enthalpy_and_boundaries<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        diffusivity_m2_s: &[f64],
        boundaries: &[SpeciesDiffusionBoundary<N>],
        duration_s: f64,
        face_partial_enthalpy: impl Fn([usize; 2], usize) -> Result<[f64; N], ChtError>,
        boundary_partial_enthalpy: impl Fn(usize) -> Result<[f64; N], ChtError>,
        controls: SpeciesDiffusionControls<N>,
        maximum_enthalpy_defect_j: f64,
        gate: &CancelGate,
    ) -> Result<BoundarySpeciesEnthalpyDiffusionStep<N>, ChtError> {
        self.diffuse_species_with_enthalpy_and_boundaries_accounted(
            domain,
            species_mass_kg,
            specific_enthalpy_j_kg,
            diffusivity_m2_s,
            boundaries,
            duration_s,
            face_partial_enthalpy,
            boundary_partial_enthalpy,
            controls,
            maximum_enthalpy_defect_j,
            gate,
        )
        .map(|(step, _)| step)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn diffuse_species_with_enthalpy_and_boundaries_accounted<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        diffusivity_m2_s: &[f64],
        boundaries: &[SpeciesDiffusionBoundary<N>],
        duration_s: f64,
        face_partial_enthalpy: impl Fn([usize; 2], usize) -> Result<[f64; N], ChtError>,
        boundary_partial_enthalpy: impl Fn(usize) -> Result<[f64; N], ChtError>,
        controls: SpeciesDiffusionControls<N>,
        maximum_enthalpy_defect_j: f64,
        gate: &CancelGate,
    ) -> Result<
        (
            BoundarySpeciesEnthalpyDiffusionStep<N>,
            BoundaryEnthalpyAccount<N>,
        ),
        ChtError,
    > {
        validate_enthalpy(
            domain,
            specific_enthalpy_j_kg,
            maximum_enthalpy_defect_j,
            gate,
        )?;
        let thermal = BoundaryEnthalpyTransfers {
            thermal: EnthalpyTransfers {
                source: face_partial_enthalpy,
                faces: Vec::new(),
                increments: Vec::new(),
                mass_frame_defect: EnergySum::default(),
            },
            source: boundary_partial_enthalpy,
            faces: Vec::new(),
            account: BoundaryEnthalpyAccount {
                increments: Vec::new(),
                faces: Vec::new(),
                external: EnergySum::default(),
                positive: EnergySum::default(),
                negative: EnergySum::default(),
                mass_into_domain: [EnergySum::default(); N],
            },
        };
        let mut observer =
            ReservoirTransfers::new(domain, self, diffusivity_m2_s, boundaries, thermal);
        let species = self.diffuse_species_observed(
            domain,
            species_mass_kg,
            diffusivity_m2_s,
            duration_s,
            controls,
            gate,
            |_| poll(gate),
            &mut observer,
        )?;
        let boundary = observer.report(gate)?;
        let thermal_boundary = observer.observer.report(boundaries, gate)?;
        let BoundaryEnthalpyTransfers {
            thermal,
            mut account,
            ..
        } = observer.observer;
        let (step, increments) = self.finish_diffusive_enthalpy(
            domain,
            specific_enthalpy_j_kg,
            species,
            thermal.increments,
            thermal.mass_frame_defect,
            account.external,
            maximum_enthalpy_defect_j,
            gate,
        )?;
        account.increments = increments;
        Ok((
            BoundarySpeciesEnthalpyDiffusionStep {
                step,
                boundary,
                thermal_boundary,
            },
            account,
        ))
    }
}
