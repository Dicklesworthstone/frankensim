//! Source-bound enthalpy on the original surface-reservoir species transfers.
//!
//! Each declared PHYSICAL face temperature supplies one common local h_k(T)
//! for both directions of the signed barycentric counterflux. A remote gas
//! reservoir's temperature does not determine this trace across a finite
//! mass-transfer film. No donor-temperature rule, Fourier heat flux, thermal
//! film law, surface reaction or additional heat of mixing is inferred.

use super::super::{BoundarySpeciesEnthalpyDiffusionStep, SpeciesDiffusionBoundary};
use super::{
    CancelGate, ChtError, GasSpeciesEnthalpyField, GasSpeciesEnthalpySource, MaterialAdvectionStep,
    SpeciesDiffusionControls, Sum, ThermochemicalSpeciesDiffusionReport, VoxelDomain, invalid,
    poll, poll_every,
};

mod adaptive;

/// Exact native boundary order paired with explicit physical face temperatures.
///
/// The full original canonical gas bank evaluates every partial enthalpy,
/// including species absent from the adjacent cell. Temperatures and reservoir
/// compositions are held over this interval. Their source owner must resolve
/// any external thermal film/interface problem before supplying these traces.
#[derive(Clone, Copy, Debug)]
pub struct GasSpeciesEnthalpyBoundarySource<'s, 'a, const N: usize> {
    source: &'s GasSpeciesEnthalpySource<'a, N>,
    boundaries: &'s [SpeciesDiffusionBoundary<N>],
    face_temperatures_k: &'s [f64],
}

impl<'s, 'a, const N: usize> GasSpeciesEnthalpyBoundarySource<'s, 'a, N> {
    /// Bind each native boundary to its physical face temperature without
    /// reordering or changing any reservoir composition or film resistance.
    pub fn new(
        source: &'s GasSpeciesEnthalpySource<'a, N>,
        boundaries: &'s [SpeciesDiffusionBoundary<N>],
        face_temperatures_k: &'s [f64],
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        if boundaries.len() != face_temperatures_k.len() {
            return Err(invalid(
                "each gas species boundary needs its own physical face temperature",
            ));
        }
        for (index, &temperature) in face_temperatures_k.iter().enumerate() {
            poll_every(gate, index)?;
            if !temperature.is_finite()
                || temperature < source.controls.minimum_temperature_k
                || temperature > source.controls.maximum_temperature_k
            {
                return Err(invalid(
                    "gas boundary temperature leaves the admitted original caloric branch",
                ));
            }
        }
        poll(gate)?;
        Ok(Self {
            source,
            boundaries,
            face_temperatures_k,
        })
    }

    /// Original cell-state owner, source bank and transport binding.
    pub fn source(&self) -> &'s GasSpeciesEnthalpySource<'a, N> {
        self.source
    }
    /// Native boundary descriptors in the exact bound order.
    pub fn boundaries(&self) -> &'s [SpeciesDiffusionBoundary<N>] {
        self.boundaries
    }
    /// Physical interface traces in the same order; these are not remote donor temperatures.
    pub fn face_temperatures_k(&self) -> &'s [f64] {
        self.face_temperatures_k
    }

    /// Refresh the SAME source calorics/transport from a complete cell candidate.
    pub fn prepare(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<GasSpeciesEnthalpyField<'a, N>, ChtError> {
        self.source.prepare(
            domain,
            carrier,
            species_mass_kg,
            specific_enthalpy_j_kg,
            gate,
        )
    }

    /// Retain the original interior-face interpolation and partial enthalpy law.
    pub fn face_partials(
        &self,
        field: &GasSpeciesEnthalpyField<'a, N>,
        domain: &VoxelDomain,
        cells: [usize; 2],
        axis: usize,
        gate: &CancelGate,
    ) -> Result<[f64; N], ChtError> {
        self.source.face_partials(field, domain, cells, axis, gate)
    }

    /// Evaluate all source partial SPECIFIC enthalpies at one declared face [J/kg].
    ///
    /// The same local values multiply positive and negative constituent
    /// transfers. Pressure-flow work already contained in h is not added again.
    /// Composition/film/face topology admission remains with the native
    /// reservoir mass owner; this query also refuses an invalid face or source
    /// binding before indexing the prepared field.
    pub fn boundary_partials(
        &self,
        field: &GasSpeciesEnthalpyField<'a, N>,
        domain: &VoxelDomain,
        boundary_index: usize,
        gate: &CancelGate,
    ) -> Result<[f64; N], ChtError> {
        poll(gate)?;
        let boundary = self
            .boundaries
            .get(boundary_index)
            .ok_or_else(|| invalid("gas caloric boundary index is outside the bound list"))?;
        if !std::ptr::eq(field.source, self.source.source)
            || !std::ptr::eq(field.context, self.source.context)
            || field.species != self.source.species
            || field.temperatures_k.len() != domain.cell_count()
            || boundary.cell >= domain.cell_count()
            || boundary.face >= 6
        {
            return Err(invalid(
                "gas boundary query does not match its prepared source/domain",
            ));
        }
        if !domain.is_fluid(boundary.cell)
            || domain
                .neighbor(boundary.cell, boundary.face)
                .is_some_and(|neighbor| domain.is_fluid(neighbor))
            || field.transport(boundary.cell).is_none()
        {
            return Err(invalid(
                "gas caloric reservoir must lie on an exterior or fluid-solid face",
            ));
        }
        self.source
            .partials_at_temperature(self.face_temperatures_k[boundary_index], gate)
    }

    /// Owned prepared arrays only; the bound temperatures/reservoirs are borrowed.
    pub fn maximum_prepared_bytes(&self, cells: usize) -> Result<usize, ChtError> {
        self.source.maximum_prepared_bytes(cells)
    }
}

/// One complete open species/enthalpy candidate and refreshed source gas state.
///
/// Signed external reservoir energy is already included in step.step's h.
/// The original boundary mass/energy accounts and cell source-inversion
/// residuals remain separate; no extra thermal deposit is required.
#[derive(Clone, Debug)]
pub struct BoundaryThermochemicalSpeciesDiffusionStep<'a, const N: usize> {
    pub step: BoundarySpeciesEnthalpyDiffusionStep<N>,
    pub state: GasSpeciesEnthalpyField<'a, N>,
    pub report: ThermochemicalSpeciesDiffusionReport,
}

impl MaterialAdvectionStep {
    /// Advance interior and authored surface diffusion with original gas calorics.
    ///
    /// Initial cell D/interior h_k and explicit boundary-face h_k are frozen
    /// for this interval. Actual boundary map transfers supply external energy
    /// once through the native open enthalpy owner. The endpoint inverts both
    /// resulting inventories together and refreshes original gas transport.
    /// Source density remains diagnostic for subsequent pressure coupling; the
    /// accepted carrier and external reservoir storage are not changed.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_thermochemical_species_with_boundaries<'s, 'a, const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        duration_s: f64,
        source: &GasSpeciesEnthalpyBoundarySource<'s, 'a, N>,
        mut controls: SpeciesDiffusionControls<N>,
        maximum_enthalpy_defect_j: f64,
        gate: &CancelGate,
    ) -> Result<BoundaryThermochemicalSpeciesDiffusionStep<'a, N>, ChtError> {
        poll(gate)?;
        if !maximum_enthalpy_defect_j.is_finite() || maximum_enthalpy_defect_j < 0.0 {
            return Err(invalid(
                "finite nonnegative gas boundary enthalpy defect budget required",
            ));
        }
        let prepared_bytes = source
            .maximum_prepared_bytes(domain.cell_count())?
            .checked_add(core::mem::size_of::<ThermochemicalSpeciesDiffusionReport>())
            .ok_or_else(|| invalid("gas boundary preparation storage extent overflows"))?;
        controls.maximum_storage_bytes = controls
            .maximum_storage_bytes
            .checked_sub(prepared_bytes)
            .ok_or_else(|| {
                invalid("gas boundary preparation exceeds the diffusion storage budget")
            })?;
        let initial =
            source.prepare(domain, self, species_mass_kg, specific_enthalpy_j_kg, gate)?;
        let initial_caloric_defect_j = initial.absolute_caloric_residual_j();
        let remaining = maximum_enthalpy_defect_j - initial_caloric_defect_j;
        if remaining < 0.0 {
            return Err(invalid(
                "initial gas boundary caloric closure exceeds the enthalpy budget",
            ));
        }
        let step = self.diffuse_species_with_enthalpy_and_boundaries(
            domain,
            species_mass_kg,
            specific_enthalpy_j_kg,
            initial.diffusivity_m2_s(),
            source.boundaries,
            duration_s,
            |cells, axis| source.face_partials(&initial, domain, cells, axis, gate),
            |index| source.boundary_partials(&initial, domain, index, gate),
            controls,
            remaining,
            gate,
        )?;
        drop(initial);
        let state = source.prepare(
            domain,
            self,
            &step.step.species.species_mass_kg,
            &step.step.specific_enthalpy_j_kg,
            gate,
        )?;
        let endpoint_caloric_defect_j = state.absolute_caloric_residual_j();
        let mut defect = Sum::default();
        for value in [
            step.step.report.absolute_enthalpy_defect_j,
            initial_caloric_defect_j,
            endpoint_caloric_defect_j,
        ] {
            defect.add(value);
        }
        if !defect.value.is_finite() || defect.value > maximum_enthalpy_defect_j {
            return Err(invalid(
                "gas boundary source closure and transport exceed the enthalpy budget",
            ));
        }
        poll(gate)?;
        Ok(BoundaryThermochemicalSpeciesDiffusionStep {
            step,
            state,
            report: ThermochemicalSpeciesDiffusionReport {
                initial_caloric_defect_j,
                endpoint_caloric_defect_j,
                absolute_enthalpy_defect_j: defect.value,
            },
        })
    }
}
