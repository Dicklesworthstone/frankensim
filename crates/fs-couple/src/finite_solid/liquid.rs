//! Release a relaxed accepted finite body into the original thermal liquid.
//!
//! The body supplies its actual reference mass, displacement and velocity.
//! The existing thermal timeline supplies every enthalpy on its original chart.
//! This zero-duration representation change is available only after the source
//! is fully liquid and its elastic stress/storage are negligible within explicit
//! bounds. Plastic, damage, viscoelastic and external contact histories do not
//! become disposable by exposing a lower-level mechanical snapshot.

mod admission;
mod balance;

use fs_exec::{CancelGate, Cx};
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceEnthalpyControls, FreeSurfaceEnthalpyError, MaterialFreeSurfaceEnthalpy3,
};
use fs_lbm::material::free_surface::transfer::enthalpy::{
    TetEnthalpyTransferError, TetEnthalpyTransferLimits, TetEnthalpyTransferReport,
    TetrahedralLiquidEnthalpyState,
};
use fs_lbm::material::free_surface::transfer::{
    LiquidTransferMoments, TetReferenceDensityTransferControls, TetrahedralLiquidReferenceDensity,
    TetrahedralLiquidState,
};
use fs_lbm::material::free_surface::{MaterialFreeSurfaceDomain3, MaterialFreeSurfaceScaling};
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};

use super::{FiniteSolidBody, FiniteSolidError};

/// Existing native budgets plus explicit source-state and representation checks.
#[derive(Clone, Copy, Debug)]
pub struct FiniteSolidLiquidReleaseControls {
    pub transfer: TetEnthalpyTransferLimits,
    pub thermal: FreeSurfaceEnthalpyControls,
    pub density: TetReferenceDensityTransferControls,
    /// Maximum Frobenius norm of ANY native current Cauchy stress [Pa].
    pub maximum_cauchy_stress_pa: f64,
    /// Sum of absolute original element elastic potentials [J]. This is an
    /// admitted unresolved storage residue, never a physical relaxation law.
    pub maximum_elastic_storage_j: f64,
    /// Maximum nodal norm of the accepted native M*a [N]. The source must have
    /// no retained load/support owner; an accelerating loaded endpoint refuses.
    pub maximum_inertial_force_n: f64,
    /// Original consistent-mass vs exact tetrahedral kinetic identity [J].
    pub source_kinetic_tolerance_j: f64,
    pub source_centroid_tolerance_m: f64,
    /// Original caloric temperature vs current held mechanical-source T [K].
    pub source_temperature_tolerance_k: f64,
    /// Sum of absolute original thermal-cell vs mechanical reference volumes [m3].
    pub source_reference_volume_tolerance_m3: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FiniteSolidLiquidBalanceStage {
    TetrahedralSource,
    NativeLiquid,
}

/// Signed receiving representation minus the original accepted source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidLiquidDifference {
    pub mass_kg: f64,
    pub momentum_kg_m_s: [f64; 3],
    pub angular_momentum_kg_m2_s: [f64; 3],
    pub centroid_m: [f64; 3],
    pub kinetic_j: f64,
    pub enthalpy_j: f64,
}

/// No error changes either borrowed accepted owner.
#[derive(Debug)]
pub enum FiniteSolidLiquidReleaseError {
    Invalid(&'static str),
    Solid(FiniteSolidError),
    Thermal(DepositionError),
    Phase(PhaseStateError),
    Transfer(TetEnthalpyTransferError),
    NativeThermal(FreeSurfaceEnthalpyError),
    Balance {
        stage: FiniteSolidLiquidBalanceStage,
        difference: FiniteSolidLiquidDifference,
    },
    Cancelled,
}

impl core::fmt::Display for FiniteSolidLiquidReleaseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-solid liquid release: {self:?}")
    }
}
impl std::error::Error for FiniteSolidLiquidReleaseError {}
impl From<FiniteSolidError> for FiniteSolidLiquidReleaseError {
    fn from(value: FiniteSolidError) -> Self {
        Self::Solid(value)
    }
}
impl From<DepositionError> for FiniteSolidLiquidReleaseError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl From<PhaseStateError> for FiniteSolidLiquidReleaseError {
    fn from(value: PhaseStateError) -> Self {
        Self::Phase(value)
    }
}
impl From<TetEnthalpyTransferError> for FiniteSolidLiquidReleaseError {
    fn from(value: TetEnthalpyTransferError) -> Self {
        Self::Transfer(value)
    }
}
impl From<FreeSurfaceEnthalpyError> for FiniteSolidLiquidReleaseError {
    fn from(value: FreeSurfaceEnthalpyError) -> Self {
        Self::NativeThermal(value)
    }
}

/// Actual native source inventories and independently measured receiver defects.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteSolidLiquidReleaseReport {
    pub physical_time_s: f64,
    /// Original persistent cells, in the original mechanical tetrahedron order.
    pub thermal_cells: Vec<CellId>,
    /// Consistent-mass native M/P/L/K at the supplied current grid origin.
    pub solid: LiquidTransferMoments,
    /// Original cached total reported by the accepted finite-body owner.
    pub reported_reference_mass_kg: f64,
    /// Sum of original thermal-cell m*h on the exact same retained chart [J].
    pub thermal_enthalpy_j: f64,
    pub source_cell_mass_defect_kg: f64,
    pub source_reference_volume_defect_m3: f64,
    pub maximum_source_temperature_defect_k: f64,
    pub maximum_cauchy_stress_pa: f64,
    pub maximum_inertial_force_n: f64,
    /// Signed original elastic storage removed from the representation [J].
    /// It must fit the near-zero gate and the complete projection budget.
    pub elastic_storage_j: f64,
    pub absolute_elastic_storage_j: f64,
    pub source_identity: FiniteSolidLiquidDifference,
    pub receiver_identity: FiniteSolidLiquidDifference,
    /// Original thermal cells vs actual native cells, summed over detached pieces.
    pub sum_absolute_component_enthalpy_defect_j: f64,
    /// All observed source/overlap/native mass discrepancies share this budget.
    pub mass_error_budget_used_kg: f64,
    pub transfer: TetEnthalpyTransferReport,
    kinetic_change: Dd,
    mechanical_change: Dd,
}

impl FiniteSolidLiquidReleaseReport {
    /// Liquid kinetic storage minus the actual original consistent-mass storage.
    pub const fn retained_kinetic_change_j(&self) -> Dd {
        self.kinetic_change
    }
    /// Kinetic projection change minus original elastic storage. Gravity is
    /// absent in this unloaded handoff. No part of this signed defect is heat.
    pub const fn retained_mechanical_change_j(&self) -> Dd {
        self.mechanical_change
    }
}

impl FiniteSolidBody<'_> {
    /// Return a native thermal-liquid candidate at this SAME accepted clock.
    ///
    /// This accepts the actual plain hyperelastic owner and complete original
    /// thermal timeline. It does not accept a plastic/damage/viscoelastic owner,
    /// discharge another coupled owner's history, or determine melting kinetics.
    /// Each thermal cell maps to one original tet, retains the same caloric
    /// object and reference mass/volume, and is already fully liquid. Shared
    /// vertices must have identical incident-cell h: no P0-to-P1 averaging or
    /// hidden heat redistribution is performed.
    ///
    /// Native stress, signed/absolute elastic storage and M*a are measured.
    /// Near-zero residues remain reported numerical discrepancies and share
    /// the existing mechanical projection budget. The operation introduces no
    /// heat, pressure work or impulse and advances neither accepted source.
    ///
    /// Stationary supports, retained pressure/gravity, initialization loads,
    /// thermal links/ports/future events and nonzero receiving gravity refuse.
    /// Those owners must first be discharged by their physical transition;
    /// they cannot silently disappear with this representation. Original
    /// temperatures, source contexts and pinned density claims must agree.
    ///
    /// On success the caller replaces BOTH original representations with this
    /// candidate. Continuing them as additional matter would duplicate mass.
    /// A refusal leaves the complete accepted source available for retry.
    #[allow(clippy::too_many_arguments)]
    pub fn release_relaxed_liquid_candidate<'s, 'm>(
        &self,
        thermal: &DepositionSimulation<'_>,
        cells_in_tet_order: &[CellId],
        scaling: MaterialFreeSurfaceScaling<'s>,
        domain: MaterialFreeSurfaceDomain3<'_>,
        grid_origin_m: [f64; 3],
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FiniteSolidLiquidReleaseControls,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> Result<
        (
            MaterialFreeSurfaceEnthalpy3<'s, 'm>,
            FiniteSolidLiquidReleaseReport,
        ),
        FiniteSolidLiquidReleaseError,
    > {
        check(cx, gate)?;
        admission::preflight(self, thermal, cells_in_tet_order, domain, controls)?;
        let positions = self.current_positions(cx)?;
        let displacement = self.displacement_m();
        let velocity = self.velocity_m_s();
        let source = admission::source(
            self,
            thermal,
            cells_in_tet_order,
            &displacement,
            scaling,
            law,
            caloric,
            controls,
            cx,
            gate,
        )?;
        let original = balance::native(self, &positions, &velocity, grid_origin_m, cx, gate)?;
        if original.maximum_inertial_force_n > controls.maximum_inertial_force_n {
            return Err(invalid(
                "accepted native inertia still carries an unresolved load",
            ));
        }
        let field = TetrahedralLiquidEnthalpyState {
            motion: TetrahedralLiquidState {
                positions_m: &positions,
                tetrahedra: &source.tetrahedra,
                velocities_m_s: &velocity,
            },
            specific_enthalpy_j_kg: &source.enthalpy,
        };
        let (receiver, transfer) = scaling.free_surface3_from_reference_tetrahedra_with_enthalpy(
            domain,
            grid_origin_m,
            field,
            TetrahedralLiquidReferenceDensity {
                positions_m: self.problem.nodes_m,
                corner_density_kg_m3: &source.density,
            },
            controls.density,
            law,
            caloric,
            controls.transfer,
            controls.thermal,
            cx,
            gate,
        )?;
        let p1 = transfer
            .mechanical
            .source_kinetic_in_frame(original.frame.reference_velocity_m_s)
            .map_err(TetEnthalpyTransferError::Transfer)?;
        let source_kinetic = balance::kinetic_change(original.frame, p1)?;
        let source_identity = balance::audit(
            original.moments,
            source.enthalpy_j.to_f64(),
            transfer.mechanical.source,
            transfer.source_enthalpy_j,
            source_kinetic,
            controls,
            FiniteSolidLiquidBalanceStage::TetrahedralSource,
        )?;
        let kinetic_change =
            source_kinetic + transfer.mechanical.retained_kinetic_energy_change_j();
        let receiver_identity = balance::audit(
            original.moments,
            source.enthalpy_j.to_f64(),
            transfer.mechanical.target,
            transfer.target_enthalpy_j,
            kinetic_change,
            controls,
            FiniteSolidLiquidBalanceStage::NativeLiquid,
        )?;
        let reference = transfer
            .mechanical
            .reference_density
            .as_ref()
            .ok_or_else(|| {
                invalid("native transfer lost the original reference-density receipt")
            })?;
        let mass_error_budget_used_kg = source.mass_defect_kg
            + (self.reference_mass_kg - original.moments.mass_kg).abs()
            + source_identity.mass_kg.abs()
            + reference.overlap_mass_partition_defect_kg
            + reference.sum_absolute_cell_mass_defect_kg;
        let mut projection = source_kinetic.abs() + source.absolute_elastic_j;
        let mut component_h_error = Dd::ZERO;
        let mut component_mass_error = Dd::ZERO;
        for (piece, heat) in transfer
            .mechanical
            .components
            .iter()
            .zip(&transfer.components)
        {
            check(cx, gate)?;
            let mut mass = Dd::ZERO;
            let mut enthalpy = Dd::ZERO;
            for &tet in &piece.source_tetrahedra {
                check(cx, gate)?;
                let cell = &thermal.cells()[&cells_in_tet_order[tet]];
                mass = mass + Dd::from_f64(cell.mass_kg);
                enthalpy =
                    enthalpy + Dd::from_f64(cell.mass_kg) * Dd::from_f64(cell.specific_h_j_kg);
            }
            component_mass_error =
                component_mass_error + (Dd::from_f64(piece.target.mass_kg) - mass).abs();
            component_h_error =
                component_h_error + (Dd::from_f64(heat.target_enthalpy_j) - enthalpy).abs();
            projection = projection + Dd::from_f64(piece.kinetic_energy_change_j).abs();
        }
        let mechanical_change = kinetic_change - source.elastic_j;
        if !mass_error_budget_used_kg.is_finite()
            || mass_error_budget_used_kg > controls.transfer.mechanical.mass_tolerance_kg
            || !finite_dd(projection)
            || !finite_dd(mechanical_change)
            || !finite_dd(component_h_error)
            || !finite_dd(component_mass_error)
            || Dd::from_f64(
                controls
                    .transfer
                    .mechanical
                    .maximum_projection_energy_change_j,
            )
            .lt(projection)
            || Dd::from_f64(controls.transfer.enthalpy_tolerance_j).lt(component_h_error)
            || Dd::from_f64(controls.transfer.mechanical.mass_tolerance_kg).lt(component_mass_error)
        {
            return Err(invalid(
                "original solid/thermal and native liquid exceed their shared conservation or mechanical projection budgets",
            ));
        }
        let receiver = receiver.with_time_origin(self.state.t)?;
        check(cx, gate)?;
        Ok((
            receiver,
            FiniteSolidLiquidReleaseReport {
                physical_time_s: self.state.t,
                thermal_cells: cells_in_tet_order.to_vec(),
                solid: original.moments,
                reported_reference_mass_kg: self.reference_mass_kg,
                thermal_enthalpy_j: source.enthalpy_j.to_f64(),
                source_cell_mass_defect_kg: source.mass_defect_kg,
                source_reference_volume_defect_m3: source.volume_defect_m3,
                maximum_source_temperature_defect_k: source.maximum_temperature_defect_k,
                maximum_cauchy_stress_pa: source.maximum_cauchy_stress_pa,
                maximum_inertial_force_n: original.maximum_inertial_force_n,
                elastic_storage_j: source.elastic_j.to_f64(),
                absolute_elastic_storage_j: source.absolute_elastic_j.to_f64(),
                source_identity,
                receiver_identity,
                sum_absolute_component_enthalpy_defect_j: component_h_error.to_f64(),
                mass_error_budget_used_kg,
                transfer,
                kinetic_change,
                mechanical_change,
            },
        ))
    }
}

fn invalid(message: &'static str) -> FiniteSolidLiquidReleaseError {
    FiniteSolidLiquidReleaseError::Invalid(message)
}
fn check(cx: &Cx<'_>, gate: &CancelGate) -> Result<(), FiniteSolidLiquidReleaseError> {
    if gate.is_requested() || cx.checkpoint().is_err() {
        Err(FiniteSolidLiquidReleaseError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
