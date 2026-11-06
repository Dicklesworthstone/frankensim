//! Held prescribed adhesive tractions on the complete original modal network.
//!
//! Each polar cell retains its original source force and first spatial moments.
//! A supplied affine modal trace and its transpose give exactly conjugate
//! displacement and generalized load, including a zero-resultant couple.
//! Native mechanics, material memory, radiation and losses use the unchanged
//! staged owner; independent spatial work is admitted before any publication.
//!
//! These are externally prescribed, held pressure fields. They neither solve
//! implicit adhesion/compliance nor add a source potential-energy account.
//! Modal traces are an authored receiving approximation, not a new pressure law
//! or a claim that arbitrary mode shapes are affine. A field and its trace use
//! the SAME fixed contact centre and local axes for the whole sample. Original
//! network external work already includes this load: never debit it again.

use super::{CoupledModalSystem, ModalAcousticTimeModel, ModalCouplingError, finite, poll};
use fs_contact::normal_patch::traction::{
    AdhesiveTractionProfile,
    sectors::{
        MAX_POLAR_TRACTION_SECTORS, MAX_RADIAL_MOMENT_EVALUATIONS, PolarRingScratch,
        PolarTractionTransfer, PolarTransferBudget, SectorAdhesiveTraction, SectorTractionError,
    },
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Hard bound on stored cell-by-mode triples, independent of native mode caps.
pub const MAX_ADHESIVE_TRACE_TERMS: usize = 1_048_576;

/// One retained mode's normal material-displacement trace on one polar cell.
/// For mass-normalized q [m sqrt(kg)], u_z=(value+gradient·[x,y])*q.
/// A two-body trace uses the signed DIFFERENCE of their material motions,
/// giving opposite efforts through the same transpose. Moving a quadrature
/// location does not contribute material velocity.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AffineNormalModalTrace {
    /// Value at the original profile centre [1/sqrt(kg)].
    pub value: f64,
    /// Spatial derivatives in the profile's local x,y axes [1/(m sqrt(kg))].
    pub gradient_per_m: [f64; 2],
}

/// Complete polar partition and immutable original-basis receiving traces.
#[derive(Clone, Debug, PartialEq)]
pub struct PolarModalTrace {
    /// Strictly increasing radii, starting at zero; cover every supplied profile.
    pub radii_m: Vec<f64>,
    /// Strictly increasing angular edges, exactly 0..2π.
    pub angles_rad: Vec<f64>,
    /// Cell-major, then original component/mode order. Cells are radial-major,
    /// angular-minor. Include exact zeros for all unattached retained modes.
    pub coefficients: Vec<AffineNormalModalTrace>,
}

/// Cold storage/work admission and contact-sized accuracy controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveTractionModalConfig {
    pub maximum_trace_terms: usize,
    pub transfer: PolarTransferBudget,
    /// Projection/work roundoff gate; never scaled by unrelated body storage
    /// or a large background external force.
    pub work_absolute_tolerance_j: f64,
    pub work_relative_tolerance: f64,
    /// Gate on the summed moment-quadrature work ERROR ESTIMATE, not a proof
    /// bound. Modal-trace approximation and material uncertainty are excluded.
    pub maximum_quadrature_work_estimate_j: f64,
}

#[derive(Debug)]
pub enum AdhesiveTractionModalError {
    Invalid(&'static str),
    Native(ModalCouplingError),
    Spatial(SectorTractionError),
    Work { residual_j: f64, tolerance_j: f64 },
    QuadratureWork { estimate_j: f64, maximum_j: f64 },
}
impl core::fmt::Display for AdhesiveTractionModalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "prescribed adhesive modal traction refused: {self:?}")
    }
}
impl std::error::Error for AdhesiveTractionModalError {}
impl From<ModalCouplingError> for AdhesiveTractionModalError {
    fn from(value: ModalCouplingError) -> Self {
        Self::Native(value)
    }
}
impl From<SectorTractionError> for AdhesiveTractionModalError {
    fn from(value: SectorTractionError) -> Self {
        Self::Spatial(value)
    }
}

/// Last accepted load account. The native frame retains all original storage,
/// losses and total external work; these are subdivisions of that SAME work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdhesiveTractionModalFrame {
    pub sample: u64,
    /// Generated directly from the original immutable profile, never supplied
    /// by the caller alongside mutable, unbound sector records.
    pub transfer: PolarTractionTransfer,
    /// Integral of held traction against the actual affine displacement [J].
    pub spatial_traction_work_j: f64,
    /// Work of the correction actually represented in the applied native force
    /// vector, after adding the caller's other external forces [J].
    pub traction_work_on_network_j: f64,
    pub other_external_work_j: f64,
    pub work_residual_j: f64,
    pub work_tolerance_j: f64,
    pub quadrature_work_estimate_j: f64,
}

/// Fixed trace mapping with reusable caller-owned adapter scratch. Construction
/// allocates all adapter vectors; stepping allocates no additional adapter
/// storage. The original native owner retains its own existing allocations.
/// Read-only native access prevents silently changing the modal basis. Consume
/// this wrapper to change geometry/traces, retaining all accepted native state.
pub struct AdhesiveTractionModalSystem {
    network: CoupledModalSystem,
    trace: PolarModalTrace,
    config: AdhesiveTractionModalConfig,
    cells: Vec<SectorAdhesiveTraction>,
    rings: Vec<PolarRingScratch>,
    sums: Vec<Dd>,
    held: Vec<f64>,
    delta: Vec<f64>,
    last: Option<AdhesiveTractionModalFrame>,
}

impl AdhesiveTractionModalSystem {
    pub fn new(
        network: CoupledModalSystem,
        trace: PolarModalTrace,
        config: AdhesiveTractionModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, AdhesiveTractionModalError> {
        poll(Some(gate))?;
        network.require_time_domain()?;
        check_config(config)?;
        let n = network.mode_count();
        let count = trace
            .radii_m
            .len()
            .checked_sub(1)
            .and_then(|r| {
                trace
                    .angles_rad
                    .len()
                    .checked_sub(1)
                    .and_then(|a| r.checked_mul(a))
            })
            .ok_or(AdhesiveTractionModalError::Invalid("polar trace extent"))?;
        if n == 0
            || count == 0
            || count > MAX_POLAR_TRACTION_SECTORS
            || count.checked_mul(n) != Some(trace.coefficients.len())
            || trace.coefficients.len() > config.maximum_trace_terms
            || trace.radii_m[0] != 0.0
            || trace.angles_rad[0] != 0.0
            || *trace.angles_rad.last().unwrap() != core::f64::consts::TAU
        {
            return Err(AdhesiveTractionModalError::Invalid(
                "complete bounded polar trace in original modal order",
            ));
        }
        for edge in trace.radii_m.windows(2).chain(trace.angles_rad.windows(2)) {
            if !edge[0].is_finite() || !edge[1].is_finite() || edge[1] <= edge[0] {
                return Err(AdhesiveTractionModalError::Invalid(
                    "finite increasing polar edges",
                ));
            }
        }
        for tile in trace.coefficients.chunks(256) {
            poll(Some(gate))?;
            if tile
                .iter()
                .any(|v| !v.value.is_finite() || v.gradient_per_m.iter().any(|x| !x.is_finite()))
            {
                return Err(AdhesiveTractionModalError::Invalid(
                    "finite affine modal traces",
                ));
            }
        }
        let rings = vec![PolarRingScratch::default(); trace.radii_m.len() - 1];
        poll(Some(gate))?;
        Ok(Self {
            network,
            trace,
            config,
            cells: vec![SectorAdhesiveTraction::default(); count],
            rings,
            sums: vec![Dd::ZERO; n],
            held: vec![0.0; n],
            delta: vec![0.0; n],
            last: None,
        })
    }

    #[must_use]
    pub const fn network(&self) -> &CoupledModalSystem {
        &self.network
    }
    #[must_use]
    pub fn into_network(self) -> CoupledModalSystem {
        self.network
    }
    #[must_use]
    pub const fn last_frame(&self) -> Option<&AdhesiveTractionModalFrame> {
        self.last.as_ref()
    }

    /// Apply the supplied original profile, held over one native sample. Its
    /// support must fit the fixed grid; its local axes must match the traces.
    /// `external` excludes this traction and follows original modal order.
    /// Source history/potential is never advanced or credited by this adapter.
    pub fn step(
        &mut self,
        profile: &AdhesiveTractionProfile,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<&AdhesiveTractionModalFrame, AdhesiveTractionModalError> {
        let accepted = self.stage(profile, external, gate)?;
        self.network.publish_staged();
        Ok(self.last.insert(accepted))
    }

    /// Admit another owner's candidate before mechanics or this account changes.
    /// The callback returns a candidate; installing it after success must be
    /// infallible. A callback refusal or late cancellation preserves every
    /// accepted native state, clock, port and previously accepted load account.
    pub fn step_with_admission<R, E>(
        &mut self,
        profile: &AdhesiveTractionProfile,
        external: &[f64],
        gate: &CancelGate,
        admit: impl FnOnce(&AdhesiveTractionModalFrame, &[ModalAcousticTimeModel]) -> Result<R, E>,
    ) -> Result<R, E>
    where
        E: From<AdhesiveTractionModalError>,
    {
        let accepted = self.stage(profile, external, gate).map_err(E::from)?;
        let candidate = admit(&accepted, self.network.staged_components())?;
        poll(Some(gate))
            .map_err(AdhesiveTractionModalError::from)
            .map_err(E::from)?;
        self.network.publish_staged();
        self.last = Some(accepted);
        Ok(candidate)
    }

    fn stage(
        &mut self,
        profile: &AdhesiveTractionProfile,
        external: &[f64],
        gate: &CancelGate,
    ) -> Result<AdhesiveTractionModalFrame, AdhesiveTractionModalError> {
        poll(Some(gate))?;
        let n = self.network.mode_count();
        if external.len() != n || external.iter().any(|x| !x.is_finite()) {
            return Err(AdhesiveTractionModalError::Invalid(
                "finite complete original external force vector",
            ));
        }
        let transfer = profile.transfer_polar_grid(
            &self.trace.radii_m,
            &self.trace.angles_rad,
            &mut self.cells,
            &mut self.rings,
            self.config.transfer,
            gate,
        )?;
        self.sums.fill(Dd::ZERO);
        for (cell, traces) in self.cells.iter().zip(self.trace.coefficients.chunks(n)) {
            for (tile, values) in traces.chunks(256).zip(self.sums.chunks_mut(256)) {
                poll(Some(gate))?;
                for (trace, sum) in tile.iter().zip(values) {
                    *sum = *sum
                        + product(cell.normal_force_n, trace.value)
                        + product(cell.first_moment_n_m[0], trace.gradient_per_m[0])
                        + product(cell.first_moment_n_m[1], trace.gradient_per_m[1]);
                }
            }
        }
        for ((held, input), sum) in self.held.iter_mut().zip(external).zip(&self.sums) {
            *held = finite(input + finite(sum.to_f64())?)?;
        }
        self.network.stage_inner(&self.held, Some(gate))?;
        // Rebuild from the actual candidate and accepted native states, not a
        // free prediction or a work value inferred from energy closure.
        for ((change, end), start) in self
            .delta
            .iter_mut()
            .zip(self.network.candidates.iter().flat_map(|m| m.states()))
            .zip(self.network.models.iter().flat_map(|m| m.states()))
        {
            *change = finite(end.displacement_m_sqrt_kg - start.displacement_m_sqrt_kg)?;
        }
        let mut actual = Dd::ZERO;
        let mut other = Dd::ZERO;
        for ((held, input), change) in self.held.iter().zip(external).zip(&self.delta) {
            // A small force erased by addition to a large external input must
            // not receive fictional work credit. Subtract represented inputs.
            actual = actual + (Dd::from_f64(*held) - Dd::from_f64(*input)) * Dd::from_f64(*change);
            other = other + product(*input, *change);
        }
        let mut spatial = Dd::ZERO;
        let mut scale = Dd::ZERO;
        let mut estimate = Dd::ZERO;
        for (cell, traces) in self.cells.iter().zip(self.trace.coefficients.chunks(n)) {
            let mut displacement = [Dd::ZERO; 3];
            for (tile, changes) in traces.chunks(256).zip(self.delta.chunks(256)) {
                poll(Some(gate))?;
                for (trace, &change) in tile.iter().zip(changes) {
                    displacement[0] = displacement[0] + product(trace.value, change);
                    displacement[1] = displacement[1] + product(trace.gradient_per_m[0], change);
                    displacement[2] = displacement[2] + product(trace.gradient_per_m[1], change);
                }
            }
            let d = [
                finite(displacement[0].to_f64())?,
                finite(displacement[1].to_f64())?,
                finite(displacement[2].to_f64())?,
            ];
            for (effort, change) in [
                cell.normal_force_n,
                cell.first_moment_n_m[0],
                cell.first_moment_n_m[1],
            ]
            .into_iter()
            .zip(d)
            {
                let work = product(effort, change);
                spatial = spatial + work;
                scale = scale + Dd::from_f64(finite(work.to_f64())?.abs());
            }
            estimate =
                estimate + product(cell.moment_error_estimate_n_m, finite(d[1].hypot(d[2]))?);
        }
        let spatial_work = finite(spatial.to_f64())?;
        let actual_work = finite(actual.to_f64())?;
        let work_residual = finite((actual - spatial).to_f64())?;
        let tolerance = finite(
            self.config.work_absolute_tolerance_j
                + self.config.work_relative_tolerance
                    * finite(scale.to_f64())?.max(actual_work.abs()),
        )?;
        if work_residual.abs() > tolerance {
            return Err(AdhesiveTractionModalError::Work {
                residual_j: work_residual,
                tolerance_j: tolerance,
            });
        }
        let quadrature_work = finite(estimate.to_f64())?;
        if quadrature_work > self.config.maximum_quadrature_work_estimate_j {
            return Err(AdhesiveTractionModalError::QuadratureWork {
                estimate_j: quadrature_work,
                maximum_j: self.config.maximum_quadrature_work_estimate_j,
            });
        }
        let native = &self.network.staged_frame;
        let native_residual =
            finite((Dd::from_f64(native.external_work_j) - actual - other).to_f64())?;
        if native_residual.abs() > native.energy_tolerance_j {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: native_residual,
                tolerance_j: native.energy_tolerance_j,
            }
            .into());
        }
        let accepted = AdhesiveTractionModalFrame {
            sample: native.sample,
            transfer,
            spatial_traction_work_j: spatial_work,
            traction_work_on_network_j: actual_work,
            other_external_work_j: finite(other.to_f64())?,
            work_residual_j: work_residual,
            work_tolerance_j: tolerance,
            quadrature_work_estimate_j: quadrature_work,
        };
        self.network.validate_staged_material_points(Some(gate))?;
        poll(Some(gate))?;
        Ok(accepted)
    }
}

fn product(a: f64, b: f64) -> Dd {
    Dd::from_f64(a) * Dd::from_f64(b)
}

fn check_config(c: AdhesiveTractionModalConfig) -> Result<(), AdhesiveTractionModalError> {
    let radial = c.transfer.radial_moment;
    if c.maximum_trace_terms == 0
        || c.maximum_trace_terms > MAX_ADHESIVE_TRACE_TERMS
        || radial.maximum_evaluations > MAX_RADIAL_MOMENT_EVALUATIONS
        || [
            c.work_absolute_tolerance_j,
            c.maximum_quadrature_work_estimate_j,
            c.transfer.force_absolute_tolerance_n,
            c.transfer.moment_absolute_tolerance_n_m,
            radial.absolute_tolerance_n_m,
            radial.relative_tolerance,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        || !c.work_relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&c.work_relative_tolerance)
    {
        return Err(AdhesiveTractionModalError::Invalid(
            "bounded trace storage and finite nonnegative accuracy controls",
        ));
    }
    Ok(())
}
