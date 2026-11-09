//! Explicit fixed-composition surface reservoirs on the original common-D maps.
//! The reservoir owns no evolved storage. Its signed species exchange is reported
//! separately from carrier advection; it must not be applied a second time.

use super::super::fractions;
use super::{
    AdaptiveSpeciesDiffusionControls, AdaptiveSpeciesDiffusionStep, CancelGate, ChtError,
    DiffusionObserver, Inventory, MaterialAdvectionStep, SpeciesDiffusionControls,
    SpeciesDiffusionStep, Sum, VoxelDomain, det, nonnegative, poll, poll_every, positive, refuse,
};

/// One explicitly authored surface reservoir, held fixed for the interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeciesDiffusionBoundary<const N: usize> {
    /// Fluid cell adjoining the reservoir. Solid-cell storage is not evolved.
    pub cell: usize,
    /// x-, x+, y-, y+, z-, z+; must face exterior or a solid, never another fluid.
    pub face: usize,
    /// Nonnegative unit-sum mass fractions, with the SAME constituent ordering.
    /// No normalization, solvent reconstruction or partition coefficient is used.
    pub reservoir_mass_fractions: [f64; N],
    /// Common surface-film resistance [m² s/kg], finite and nonnegative.
    /// Zero prescribes composition at the actual face, not the cell centre.
    /// The law is j_k,out = (Y_k,surface - Y_k,reservoir) / resistance.
    /// For a transfer velocity k_m, supply 1/(rho_film*k_m) using its own density.
    pub film_resistance_m2_s_kg: f64,
}

/// Integrated actual transfers at one declared face, in canonical input order.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesBoundaryTransfer<const N: usize> {
    pub cell: usize,
    pub face: usize,
    /// Half-cell plus film conductance [kg/s]; zero when cell D is zero.
    pub conductance_kg_s: f64,
    /// Signed mass received by the fluid cell [kg]. Independently accumulated,
    /// not obtained by subtracting potentially large gross totals.
    pub mass_into_cell_kg: [f64; N],
    /// Gross inward and outward constituent masses [kg]; both may be positive.
    pub incoming_mass_kg: [f64; N],
    pub outgoing_mass_kg: [f64; N],
}

/// External transfer account used in the original constituent-mass admission.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesBoundaryDiffusionReport<const N: usize> {
    pub faces: Vec<SpeciesBoundaryTransfer<N>>,
    /// Signed original transfer sum [kg], positive into the fluid domain.
    pub mass_into_domain_kg: [f64; N],
    /// Sum over every boundary map of abs(sum_k dm_k) [kg].
    /// Charged to the original carrier-partition budget, never hidden by
    /// cancellation between species, faces or successive half sweeps.
    pub absolute_carrier_transfer_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BoundarySpeciesDiffusionStep<const N: usize> {
    /// Fine physical inventories and mass-balance defects including reservoirs.
    pub step: SpeciesDiffusionStep<N>,
    pub boundary: SpeciesBoundaryDiffusionReport<N>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveBoundarySpeciesDiffusionStep<const N: usize> {
    /// Its maximum_error_ratio includes both concentrations and face transfers.
    pub step: AdaptiveSpeciesDiffusionStep<N>,
    pub boundary: SpeciesBoundaryDiffusionReport<N>,
    /// Per species, sum over faces of the largest signed/inward/outward
    /// coarse/fine transfer difference divided by three [kg].
    pub estimated_boundary_transfer_error_kg: [f64; N],
}

struct BoundaryMap<const N: usize> {
    input_index: usize,
    conductance: f64,
    relaxation: f64,
    signed: [Sum; N],
    incoming: [Sum; N],
    outgoing: [Sum; N],
}

/// Optional consumer of the SAME original reservoir maps. The unit consumer
/// compiles to the existing mass-only path without thermal allocation or work.
pub(super) trait ReservoirObserver<const N: usize>: DiffusionObserver<N> {
    fn reservoir_storage_bytes(&self, _: usize, _: usize) -> Option<usize> {
        Some(0)
    }
    fn prepare_reservoir_storage(&mut self, _: usize, _: usize) {}
    fn prepare_reservoir(
        &mut self,
        _: usize,
        _: &SpeciesDiffusionBoundary<N>,
        _: &CancelGate,
    ) -> Result<(), ChtError> {
        Ok(())
    }
    fn reservoir_transferred(&mut self, _: usize, _: [f64; N]) -> Result<(), ChtError> {
        Ok(())
    }
}

impl<const N: usize> ReservoirObserver<N> for () {}

pub(super) struct ReservoirTransfers<'a, const N: usize, O = ()> {
    domain: &'a VoxelDomain,
    carrier: &'a MaterialAdvectionStep,
    diffusivity: &'a [f64],
    boundaries: &'a [SpeciesDiffusionBoundary<N>],
    active_faces: usize,
    maps: Vec<BoundaryMap<N>>,
    signed: [Sum; N],
    carrier_defect: Sum,
    maximum_relaxation: f64,
    pub(super) observer: O,
}

impl<const N: usize, O: ReservoirObserver<N>> DiffusionObserver<N>
    for ReservoirTransfers<'_, N, O>
{
    fn boundary_faces(&mut self, gate: &CancelGate) -> Result<usize, ChtError> {
        let mut previous = None;
        for (index, boundary) in self.boundaries.iter().enumerate() {
            poll_every(gate, index)?;
            let key = (boundary.cell, boundary.face);
            if boundary.cell >= self.domain.cell_count()
                || boundary.face >= 6
                || !self.domain.is_fluid(boundary.cell)
                || self
                    .domain
                    .neighbor(boundary.cell, boundary.face)
                    .is_some_and(|n| self.domain.is_fluid(n))
            {
                return Err(refuse(
                    "reservoir requires an exterior or fluid-solid face of a fluid cell",
                ));
            }
            if previous.is_some_and(|last| last >= key) {
                return Err(refuse(
                    "reservoir faces must be unique and ordered by (cell, face)",
                ));
            }
            previous = Some(key);
            fractions(&boundary.reservoir_mass_fractions)?;
            nonnegative(boundary.film_resistance_m2_s_kg)?;
            if self.diffusivity[boundary.cell] > 0.0 {
                self.active_faces = self
                    .active_faces
                    .checked_add(1)
                    .ok_or_else(|| refuse("reservoir face extent overflows"))?;
            }
        }
        Ok(self.active_faces)
    }

    fn storage_bytes(&self, cells: usize, faces: usize) -> Option<usize> {
        self.active_faces
            .checked_mul(core::mem::size_of::<BoundaryMap<N>>())?
            .checked_add(
                self.boundaries
                    .len()
                    .checked_mul(core::mem::size_of::<SpeciesBoundaryTransfer<N>>())?,
            )?
            .checked_add(core::mem::size_of::<Self>())?
            .checked_add(core::mem::size_of::<BoundarySpeciesDiffusionStep<N>>())?
            .checked_add(4 * core::mem::size_of::<[f64; N]>())?
            .checked_add(self.observer.storage_bytes(cells, faces)?)?
            .checked_add(
                self.observer
                    .reservoir_storage_bytes(self.boundaries.len(), self.active_faces)?,
            )
    }

    fn prepare_storage(&mut self, cells: usize, faces: usize) {
        self.maps.reserve_exact(self.active_faces);
        self.observer.prepare_storage(cells, faces);
        self.observer
            .prepare_reservoir_storage(self.boundaries.len(), self.active_faces);
    }

    fn prepare(&mut self, cells: [usize; 2], axis: usize) -> Result<(), ChtError> {
        self.observer.prepare(cells, axis)
    }

    fn transferred(&mut self, face: usize, transfers: [f64; N]) -> Result<(), ChtError> {
        self.observer.transferred(face, transfers)
    }

    fn prepare_boundaries(&mut self, half_dt: f64, gate: &CancelGate) -> Result<(), ChtError> {
        for (input_index, boundary) in self.boundaries.iter().enumerate() {
            poll_every(gate, input_index)?;
            let cell = boundary.cell;
            if self.diffusivity[cell] == 0.0 {
                continue;
            }
            let axis = boundary.face / 2;
            let density = positive(self.carrier.mass_kg[cell] / self.domain.volume(cell))?;
            let half = positive(
                (0.5 * self.domain.widths(cell)[axis])
                    / positive(density * self.diffusivity[cell])?,
            )?;
            let film = boundary.film_resistance_m2_s_kg;
            // Scale the series sum to avoid overflowing two finite resistances.
            let largest = half.max(film);
            let conductance = positive(
                (self.domain.face_area(cell, axis) / largest) / (1.0 + half.min(film) / largest),
            )?;
            let mass = self.carrier.mass_kg[cell];
            let mut exponent = (conductance / mass) * half_dt;
            if !exponent.is_normal() {
                exponent = det::exp(det::ln(conductance) - det::ln(mass) + det::ln(half_dt));
            }
            let relaxation = positive(-det::expm1(-exponent))?;
            if relaxation > 1.0 {
                return Err(refuse(
                    "reservoir relaxation exceeds its positive convex map",
                ));
            }
            self.maximum_relaxation = self.maximum_relaxation.max(relaxation);
            self.observer
                .prepare_reservoir(input_index, boundary, gate)?;
            self.maps.push(BoundaryMap {
                input_index,
                conductance,
                relaxation,
                signed: [Sum::default(); N],
                incoming: [Sum::default(); N],
                outgoing: [Sum::default(); N],
            });
        }
        Ok(())
    }

    fn boundary_half_sweep(
        &mut self,
        inventories: &mut [Inventory<N>],
        reverse: bool,
        gate: &CancelGate,
    ) -> Result<(), ChtError> {
        for position in 0..self.maps.len() {
            poll_every(gate, position)?;
            let index = if reverse {
                self.maps.len() - 1 - position
            } else {
                position
            };
            let map = &mut self.maps[index];
            let boundary = &self.boundaries[map.input_index];
            let inventory = &mut inventories[boundary.cell];
            let transfers = reservoir_species_transfers(
                inventory.mass.map(|s| s.value),
                self.carrier.mass_kg[boundary.cell],
                boundary.reservoir_mass_fractions,
                map.relaxation,
            )
            .map_err(refuse)?;
            let mut carrier_transfer = Sum::default();
            for (k, &dm) in transfers.iter().enumerate() {
                inventory.mass[k].add(dm);
                inventory.increment[k].add(dm);
                nonnegative(inventory.mass[k].value)?;
                map.signed[k].add(dm);
                self.signed[k].add(dm);
                carrier_transfer.add(dm);
                if dm >= 0.0 {
                    map.incoming[k].add(dm);
                } else {
                    map.outgoing[k].add(-dm);
                }
                for sum in [
                    inventory.mass[k],
                    inventory.increment[k],
                    map.signed[k],
                    map.incoming[k],
                    map.outgoing[k],
                    self.signed[k],
                ] {
                    if !sum.value.is_finite() || !sum.compensation.is_finite() {
                        return Err(refuse(
                            "reservoir transfer or inventory arithmetic overflowed",
                        ));
                    }
                }
            }
            self.carrier_defect.add(carrier_transfer.value.abs());
            nonnegative(self.carrier_defect.value)?;
            self.observer.reservoir_transferred(index, transfers)?;
        }
        Ok(())
    }

    fn boundary_mass_transfer(&self) -> [f64; N] {
        self.signed.map(|s| s.value)
    }
    fn boundary_carrier_defect(&self) -> f64 {
        self.carrier_defect.value
    }
    fn maximum_boundary_relaxation(&self) -> f64 {
        self.maximum_relaxation
    }
}

impl<'a, const N: usize, O: ReservoirObserver<N>> ReservoirTransfers<'a, N, O> {
    pub(super) fn new(
        domain: &'a VoxelDomain,
        carrier: &'a MaterialAdvectionStep,
        diffusivity: &'a [f64],
        boundaries: &'a [SpeciesDiffusionBoundary<N>],
        observer: O,
    ) -> Self {
        Self {
            domain,
            carrier,
            diffusivity,
            boundaries,
            observer,
            active_faces: 0,
            maps: Vec::new(),
            signed: [Sum::default(); N],
            carrier_defect: Sum::default(),
            maximum_relaxation: 0.0,
        }
    }

    pub(super) fn report(
        &self,
        gate: &CancelGate,
    ) -> Result<SpeciesBoundaryDiffusionReport<N>, ChtError> {
        let mut faces = Vec::with_capacity(self.boundaries.len());
        for (index, boundary) in self.boundaries.iter().enumerate() {
            poll_every(gate, index)?;
            faces.push(SpeciesBoundaryTransfer {
                cell: boundary.cell,
                face: boundary.face,
                conductance_kg_s: 0.0,
                mass_into_cell_kg: [0.0; N],
                incoming_mass_kg: [0.0; N],
                outgoing_mass_kg: [0.0; N],
            });
        }
        for (index, map) in self.maps.iter().enumerate() {
            poll_every(gate, index)?;
            let output = &mut faces[map.input_index];
            output.conductance_kg_s = map.conductance;
            output.mass_into_cell_kg = map.signed.map(|s| s.value);
            output.incoming_mass_kg = map.incoming.map(|s| s.value);
            output.outgoing_mass_kg = map.outgoing.map(|s| s.value);
        }
        poll(gate)?;
        Ok(SpeciesBoundaryDiffusionReport {
            faces,
            mass_into_domain_kg: self.boundary_mass_transfer(),
            absolute_carrier_transfer_kg: self.carrier_defect.value,
        })
    }
}

/// Exact affine cell-to-fixed-reservoir map, evaluated in inventory units.
/// The binary path retains complementary trace transport even when the majority
/// fractions both round to one. Inputs are admitted by the original owner.
pub(super) fn reservoir_species_transfers<const N: usize>(
    current: [f64; N],
    carrier: f64,
    reservoir: [f64; N],
    relaxation: f64,
) -> Result<[f64; N], &'static str> {
    let mut target = [0.0; N];
    for k in 0..N {
        target[k] = carrier * reservoir[k];
        if !target[k].is_finite() || (reservoir[k] > 0.0 && target[k] == 0.0) {
            return Err("positive reservoir constituent inventory is unrepresentable");
        }
    }
    let mut transfer = [0.0; N];
    if N == 2 {
        let component =
            usize::from(current[0] / carrier + reservoir[0] > current[1] / carrier + reservoir[1]);
        let drive = target[component] - current[component];
        let dm = relaxation * drive;
        if !dm.is_finite() || (drive != 0.0 && dm == 0.0) {
            return Err("nonzero reservoir species transfer is unrepresentable");
        }
        transfer[component] = dm;
        transfer[1 - component] = -dm;
    } else {
        for k in 0..N {
            let drive = target[k] - current[k];
            transfer[k] = relaxation * drive;
            if !transfer[k].is_finite() || (drive != 0.0 && transfer[k] == 0.0) {
                return Err("nonzero reservoir species transfer is unrepresentable");
            }
        }
    }
    Ok(transfer)
}

impl MaterialAdvectionStep {
    /// Diffuse through explicitly authored surface reservoirs and interior faces.
    ///
    /// Boundaries must be in strictly increasing (cell, face) order. Unlisted
    /// exterior/fluid-solid faces stay impermeable. A listed face replaces that
    /// diffusion condition only; it does not author bulk inflow, change carrier
    /// mass, consume solid species, or modify an advective inlet composition.
    ///
    /// G = A / (w/(2*rho*D) + R_film) [kg/s]. Every constituent shares that G.
    /// Zero D blocks exchange even with zero film resistance. Each reservoir
    /// half map solves M*dY/dt = G*(Y_reservoir-Y) exactly. Canonical forward
    /// boundary/interior and reverse interior/boundary half sweeps retain the
    /// original symmetric second-order frozen-coefficient composition.
    ///
    /// Actual signed transfers close each species' mass balance. Their summed
    /// carrier-flux roundoff is included in maximum_partition_defect_kg. Work
    /// counts include both interior and reservoir maps; storage includes prepared
    /// boundary coefficients, independent transfer sums and all face readbacks.
    /// This is an isothermal common-D reservoir, not a selective membrane law,
    /// finite solid/reservoir storage, equilibrium partition model or heat owner.
    /// Failure/cancellation publishes no partial result and changes no input.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_species_with_boundaries<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        diffusivity_m2_s: &[f64],
        boundaries: &[SpeciesDiffusionBoundary<N>],
        duration_s: f64,
        controls: SpeciesDiffusionControls<N>,
        gate: &CancelGate,
    ) -> Result<BoundarySpeciesDiffusionStep<N>, ChtError> {
        let mut observer = ReservoirTransfers::new(domain, self, diffusivity_m2_s, boundaries, ());
        let step = self.diffuse_species_observed(
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
        Ok(BoundarySpeciesDiffusionStep { step, boundary })
    }

    /// Refine the SAME complete reservoir/interior interval from original masses.
    ///
    /// In addition to the original per-cell concentration criterion, each
    /// species has a positive absolute integrated-transfer tolerance [kg].
    /// For every face compare signed, gross inward and gross outward transfers;
    /// sum their largest absolute coarse/fine differences over faces and divide
    /// by three. Thus opposing boundaries or reversals cannot hide unresolved
    /// exchange behind a small net domain mass change. Both estimates are
    /// asymptotic second-order estimates, not certified spatial/model bounds.
    /// No extrapolation modifies the accepted fine inventories or transfers.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_species_with_boundaries_adaptive<const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        diffusivity_m2_s: &[f64],
        boundaries: &[SpeciesDiffusionBoundary<N>],
        duration_s: f64,
        controls: AdaptiveSpeciesDiffusionControls<N>,
        absolute_boundary_transfer_tolerance_kg: [f64; N],
        gate: &CancelGate,
    ) -> Result<AdaptiveBoundarySpeciesDiffusionStep<N>, ChtError> {
        poll(gate)?;
        for tolerance in absolute_boundary_transfer_tolerance_kg {
            positive(tolerance)?;
        }
        let extra_retained_bytes = boundaries
            .len()
            .checked_mul(core::mem::size_of::<SpeciesBoundaryTransfer<N>>())
            .and_then(|n| n.checked_add(core::mem::size_of::<SpeciesBoundaryDiffusionReport<N>>()))
            .and_then(|n| {
                n.checked_add(core::mem::size_of::<AdaptiveBoundarySpeciesDiffusionStep<N>>())
            })
            .and_then(|n| n.checked_add(4 * core::mem::size_of::<[f64; N]>() + 512))
            .ok_or_else(|| refuse("adaptive reservoir comparison storage extent overflows"))?;
        let mut estimated_boundary_transfer_error_kg = [0.0; N];
        let (step, boundary) = self.refine_species(
            domain,
            species_mass_kg,
            controls,
            gate,
            |_| poll(gate),
            extra_retained_bytes,
            |solver| {
                let result = self.diffuse_species_with_boundaries(
                    domain,
                    species_mass_kg,
                    diffusivity_m2_s,
                    boundaries,
                    duration_s,
                    solver,
                    gate,
                )?;
                Ok((result.step, result.boundary))
            },
            |coarse, fine| {
                let mut errors = [Sum::default(); N];
                for (index, (a, b)) in coarse.faces.iter().zip(&fine.faces).enumerate() {
                    poll_every(gate, index)?;
                    for k in 0..N {
                        let error = (b.mass_into_cell_kg[k] - a.mass_into_cell_kg[k])
                            .abs()
                            .max((b.incoming_mass_kg[k] - a.incoming_mass_kg[k]).abs())
                            .max((b.outgoing_mass_kg[k] - a.outgoing_mass_kg[k]).abs())
                            / 3.0;
                        errors[k].add(error);
                    }
                }
                let mut maximum_ratio = 0.0_f64;
                for k in 0..N {
                    let error = errors[k].value;
                    nonnegative(error)?;
                    estimated_boundary_transfer_error_kg[k] = error;
                    maximum_ratio =
                        maximum_ratio.max(error / absolute_boundary_transfer_tolerance_kg[k]);
                }
                Ok(maximum_ratio)
            },
        )?;
        Ok(AdaptiveBoundarySpeciesDiffusionStep {
            step,
            boundary,
            estimated_boundary_transfer_error_kg,
        })
    }
}
