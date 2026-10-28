//! The accepted complete melt enters native VOF with its original h and motion.
//!
//! On a shared cubic grid, averaging the two MAC faces onto their cell preserves
//! linear and angular momentum. Its kinetic projection deficit is bounded and
//! reported, NEVER deposited as heat. This is an explicit withdrawal of the
//! stationary mould and its boundary schedule, not moving-wall physics. New
//! permanent walls and gravity are declared by the receiving domain.

use super::*;
use crate::conjugate::FvBoundary;
use crate::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3, MaterialFreeSurfaceScaling,
    enthalpy::{
        FreeSurfaceEnthalpyControls, FreeSurfaceEnthalpyError, MaterialFreeSurfaceEnthalpy3,
        rheology::RheologySolve,
    },
    transfer::LiquidTransferMoments,
};
use fs_exec::Cx;
use fs_material::liquid::surface::LiquidCapillarity;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, transport::ResolvedPhaseTransportState,
};

/// Explicit same-time removal of stationary support and its old boundary drives.
/// Source cells must already be liquid; no solid stress or surface energy is erased.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StationaryConstraintRelease {
    /// Retire the old fixed-domain boundary schedule; use the declared VOF domain.
    Withdraw,
}

/// Absolute conservation allowances and an independently bounded representation loss.
#[derive(Clone, Copy, Debug)]
pub struct LiquidReleaseControls {
    pub mass_tolerance_kg: f64,
    pub enthalpy_tolerance_j: f64,
    pub momentum_tolerance_kg_m_s: f64,
    pub angular_momentum_tolerance_kg_m2_s: f64,
    /// Maximum native MAC kinetic energy lost by cell-centre projection, not heat.
    pub maximum_kinetic_projection_deficit_j: f64,
    /// Floating-point allowance for the native kinetic comparison, not a loss budget.
    pub kinetic_roundoff_tolerance_j: f64,
}

/// Actual original face moments and initialized native receiving-cell moments.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidReleaseReport {
    pub time_origin_s: f64,
    pub source: LiquidTransferMoments,
    pub target: LiquidTransferMoments,
    pub source_enthalpy_j: f64,
    pub target_enthalpy_j: f64,
    /// Source native K minus actual receiving K. Tiny negative roundoff is retained.
    pub kinetic_projection_deficit_j: f64,
}

/// No caller-owned thermal or native momentum state changes on refusal.
#[derive(Debug)]
pub enum LiquidReleaseError {
    Source(MaterialFlowError<DepositionError>),
    Receiver(FreeSurfaceEnthalpyError),
    Invalid(&'static str),
}
impl core::fmt::Display for LiquidReleaseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete-liquid release refused: {self:?}")
    }
}
impl std::error::Error for LiquidReleaseError {}
impl From<MaterialFlowError<DepositionError>> for LiquidReleaseError {
    fn from(error: MaterialFlowError<DepositionError>) -> Self {
        Self::Source(error)
    }
}
impl From<FreeSurfaceEnthalpyError> for LiquidReleaseError {
    fn from(error: FreeSurfaceEnthalpyError) -> Self {
        Self::Receiver(error)
    }
}

impl<R> MaterialFlowFeedback<'_, ResolvedPhaseTransportState, R, (), PhaseFlowTransport>
where
    R: Clone
        + FnMut(usize, &DepositedCell, f64) -> Result<ResolvedPhaseTransportState, MaterialError>,
{
    /// Release an accepted, completely molten fixed-grid body into native VOF.
    ///
    /// The SAME cubic cells, original registered caloric object, material card,
    /// query context and selected transport claims supply mass, h and velocity.
    /// Native retained momentum and its face-mass kinetic audit are mandatory;
    /// the applied old thermal field is never mistaken for accepted motion.
    /// No extra thermal body or coupled observer history may be discarded.
    /// The ORIGINAL receiver owns the accepted physical time and retains all
    /// existing finite-gas, wall, electrical, magnetic and heat entry points.
    ///
    /// `Withdraw` ends the stationary mould, its future schedule and external
    /// drives at this exact accepted time. No time, impulse or heat is added.
    /// The caller declares the receiving stationary walls/gravity; removed
    /// constraints are not a moving-solid calculation. Old thermal contacts end
    /// here: use the receiver's original conduction API on current wet faces.
    /// Source rheology uses `release_fully_liquid_with_rheology`; the separate
    /// phase-conductivity multiplier and buoyancy still refuse. Surface tension must be explicitly neglected:
    /// no surface-creation energy, residual elastic stress or latent heat is invented.
    ///
    /// MAC-to-cell averaging preserves P and angular P on this shared grid.
    /// Its bounded kinetic deficit remains representation error, not heat or
    /// external work. The VOF equilibrium initializer uses zero numerical pressure;
    /// it is not a pressure projection or an initially equilibrated free surface.
    /// Refine grid/time and bound the projection loss for the intended regime.
    #[allow(clippy::too_many_arguments)]
    pub fn release_fully_liquid<'s, 'm>(
        &self,
        thermal: &DepositionSimulation<'_>,
        scaling: MaterialFreeSurfaceScaling<'s>,
        domain: MaterialFreeSurfaceDomain3<'_>,
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        withdrawal: StationaryConstraintRelease,
        limits: LiquidReleaseControls,
        transport: FreeSurfaceEnthalpyControls,
        cx: &Cx<'_>,
    ) -> Result<(MaterialFreeSurfaceEnthalpy3<'s, 'm>, LiquidReleaseReport), LiquidReleaseError>
    {
        self.release_liquid(
            thermal, scaling, domain, law, caloric, withdrawal, limits, transport, None, cx,
        )
    }

    /// Carry the original phase/native shear law into the SAME moving-liquid owner.
    ///
    /// The source must already own native rheology. Its accepted physical MAC
    /// rate supplies each comparison against the target's pinned T/rate bundle;
    /// a different table, context, rate axis or unowned coefficient refuses.
    /// The original VOF constitutive solve then evaluates that SAME source at
    /// its actual force-corrected native stress, including after heat/advection.
    /// No reference shear, viscosity or temperature is frozen at the transition.
    ///
    /// This retains `release_fully_liquid`'s mass, h, momentum, kinetic-deficit
    /// and clock checks. Its equilibrium population initializer rebuilds the
    /// viscous stress on the receiving lattice; it does not transfer the MAC
    /// stress tensor. Generalized-Newtonian response has no stored elastic
    /// memory to erase, but the stress startup and time split require refinement.
    /// Phase-conductivity multipliers, buoyancy and surface energy remain named
    /// unsupported transitions, rather than disappearing from the material law.
    #[allow(clippy::too_many_arguments)]
    pub fn release_fully_liquid_with_rheology<'s, 'm>(
        &self,
        thermal: &DepositionSimulation<'_>,
        scaling: MaterialFreeSurfaceScaling<'s>,
        domain: MaterialFreeSurfaceDomain3<'_>,
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        withdrawal: StationaryConstraintRelease,
        limits: LiquidReleaseControls,
        transport: FreeSurfaceEnthalpyControls,
        rheology: RheologySolve,
        cx: &Cx<'_>,
    ) -> Result<(MaterialFreeSurfaceEnthalpy3<'s, 'm>, LiquidReleaseReport), LiquidReleaseError>
    {
        self.release_liquid(
            thermal,
            scaling,
            domain,
            law,
            caloric,
            withdrawal,
            limits,
            transport,
            Some(rheology),
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn release_liquid<'s, 'm>(
        &self,
        thermal: &DepositionSimulation<'_>,
        scaling: MaterialFreeSurfaceScaling<'s>,
        domain: MaterialFreeSurfaceDomain3<'_>,
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        withdrawal: StationaryConstraintRelease,
        limits: LiquidReleaseControls,
        transport: FreeSurfaceEnthalpyControls,
        rheology: Option<RheologySolve>,
        cx: &Cx<'_>,
    ) -> Result<(MaterialFreeSurfaceEnthalpy3<'s, 'm>, LiquidReleaseReport), LiquidReleaseError>
    {
        let _ = withdrawal;
        let refused = LiquidReleaseError::Invalid;
        self.check_map(thermal, &mut || Ok(()))?;
        self.check_unsteady_start(thermal, &mut || Ok(()))?;
        checkpoint(cx, self.gate)?;
        for value in [
            limits.mass_tolerance_kg,
            limits.enthalpy_tolerance_j,
            limits.momentum_tolerance_kg_m_s,
            limits.angular_momentum_tolerance_kg_m2_s,
            limits.maximum_kinetic_projection_deficit_j,
            limits.kinetic_roundoff_tolerance_j,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(refused("invalid release allowance"));
            }
        }
        if !self.domain.is_uniform()
            || self.domain.dims() != domain.cells
            || self.domain.dx().to_bits() != scaling.cell_size_m().to_bits()
            || self.density_kg_m3.to_bits() != scaling.material().density_kg_m3().to_bits()
            || !core::ptr::addr_eq(
                self.input.caloric as *const dyn fs_conduction::transient::deposition::ThermalLaw,
                caloric as *const EquilibriumEnthalpyPhaseCurve,
            )
            || self.endpoint.len() != thermal.cells().len()
            || scaling.material().capillarity() != LiquidCapillarity::Neglected
            || self.config.faces.iter().any(|face| {
                !matches!(face,
                FvBoundary::Wall { velocity } if *velocity == [0.0; 3])
            })
        {
            return Err(refused(
                "release needs one complete cubic liquid, stationary walls, original law and explicitly neglected capillarity",
            ));
        }
        if self.buoyancy.is_some() || self.phase_conduction.is_some() {
            return Err(refused(
                "release must retain phase-conductivity and buoyancy through explicit receiving owners",
            ));
        }
        let rate_axis = match (
            self.endpoint_shear_rate_per_s(),
            rheology,
            scaling.material().rheology_axes(),
        ) {
            (None, None, None) => None,
            (Some(rates), Some(_), Some((temperature, rate)))
                if temperature == law.temperature_axis()
                    && rates.len() == self.domain.cell_count() =>
            {
                Some(rate)
            }
            _ => {
                return Err(refused(
                    "release must retain the original native shear owner and its explicitly evolving target source",
                ));
            }
        };
        let flow = self
            .unsteady_flow()
            .ok_or(refused("release needs accepted native momentum"))?;
        let native = self
            .unsteady_heating_audit()
            .ok_or(refused("release needs the original native kinetic audit"))?;
        let mut source = LiquidTransferMoments::default();
        let mut source_h = 0.0;
        let mut initial = Vec::with_capacity(self.endpoint.len());
        let mut enthalpy = BTreeMap::new();
        let dx = self.domain.dx();
        let area = dx * dx;
        for sample in &self.endpoint {
            checkpoint(cx, self.gate)?;
            let cell = &thermal.cells()[&sample.cell];
            let phase = caloric
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(FreeSurfaceEnthalpyError::from)?;
            if phase.phase() != SolidLiquidPhase::Liquid || &phase != sample.material.phase_state()
            {
                return Err(refused(
                    "every original cell must be exactly liquid at its accepted enthalpy",
                ));
            }
            let liquid = if let Some(axis) = rate_axis {
                let rate = self
                    .endpoint_shear_rate_per_s()
                    .expect("admitted native rate")[sample.voxel];
                let point = law
                    .reference_point()
                    .clone()
                    .with(law.temperature_axis(), sample.temperature_k)
                    .and_then(|point| point.with(axis, rate))
                    .map_err(|_| {
                        refused("invalid original native temperature/shear coordinates")
                    })?;
                law.reference_state().at_query_point(law.card(), &point)
            } else {
                law.at_temperature(fs_qty::Temperature::new(sample.temperature_k))
            }
            .map_err(FreeSurfaceEnthalpyError::from)?;
            let original = sample.material.material();
            let current = liquid.material();
            if current.card_identity() != original.card_identity()
                || current.query_point() != original.query_point()
                || current.axis_quantities() != original.axis_quantities()
                || current.properties().iter().any(|property| {
                    original
                        .property(property.requirement().name())
                        .is_none_or(|old| {
                            old.requirement() != property.requirement()
                                || old.answer().receipt.selected
                                    != property.answer().receipt.selected
                                || old.value_si().to_bits() != property.value_si().to_bits()
                        })
                })
            {
                return Err(refused(
                    "release must retain every original liquid transport source and query context",
                ));
            }
            let at = self.domain.coords(sample.voxel);
            let center = self.domain.center(at[0], at[1], at[2]);
            source.mass_kg += cell.mass_kg;
            for a in 0..3 {
                source.first_mass_moment_kg_m[a] += cell.mass_kg * center[a];
            }
            for axis in 0..3 {
                for side in 0..2 {
                    let mut face = at;
                    face[axis] += side;
                    let speed = match axis {
                        0 => flow.field.flux_x(face[0], face[1], face[2]),
                        1 => flow.field.flux_y(face[0], face[1], face[2]),
                        _ => flow.field.flux_z(face[0], face[1], face[2]),
                    } / area;
                    let mut position = center;
                    position[axis] = self.domain.face_coord(axis, face[axis]);
                    let mut velocity = [0.0; 3];
                    velocity[axis] = speed;
                    accumulate_motion(&mut source, 0.5 * cell.mass_kg, position, velocity);
                }
            }
            source_h += cell.mass_kg * cell.specific_h_j_kg;
            enthalpy.insert(at, cell.specific_h_j_kg);
            initial.push(MaterialFreeSurfaceInitialCell3 {
                coordinates: at,
                fill_fraction: 1.0,
                velocity_m_s: flow.velocity_m_s[sample.voxel],
            });
        }
        if (source.kinetic_energy_j - native.final_kinetic_energy_j).abs()
            > limits.kinetic_roundoff_tolerance_j
        {
            return Err(refused(
                "source face quadrature differs from the original kinetic owner",
            ));
        }
        // The source's accepted native audit, including any supported capture,
        // remains authoritative; the independently summed faces checked it above.
        source.kinetic_energy_j = native.final_kinetic_energy_j;
        let initial_h = |at| enthalpy.get(&at).copied().unwrap_or(f64::NAN);
        let liquid = if let Some(solve) = rheology {
            scaling.free_surface3_from_cells_with_rheology(
                domain, &initial, law, caloric, transport, solve, initial_h, cx, self.gate,
            )?
        } else {
            scaling
                .free_surface3_from_cells(domain, &initial, cx)
                .map_err(FreeSurfaceEnthalpyError::from)?
                .with_enthalpy(law, caloric, transport, initial_h, self.gate)?
        }
        .with_time_origin(thermal.time_s())?;
        let mut target = LiquidTransferMoments::default();
        for cell in &initial {
            checkpoint(cx, self.gate)?;
            let [x, y, z] = cell.coordinates;
            let sample = liquid
                .flow()
                .sample_at_si(x, y, z)
                .map_err(FreeSurfaceEnthalpyError::from)?
                .ok_or(refused("native receiver lost a transferred liquid cell"))?;
            let mass = sample.fill_fraction * sample.density_kg_m3 * dx * dx * dx;
            let center = self.domain.center(x, y, z);
            target.mass_kg += mass;
            for a in 0..3 {
                target.first_mass_moment_kg_m[a] += mass * center[a];
            }
            accumulate_motion(&mut target, mass, center, sample.velocity_m_s);
        }
        let target_h = liquid.total_enthalpy_j()?;
        let deficit = source.kinetic_energy_j - target.kinetic_energy_j;
        let finite = |m: &LiquidTransferMoments| {
            [m.mass_kg, m.kinetic_energy_j]
                .into_iter()
                .chain(m.first_mass_moment_kg_m)
                .chain(m.momentum_kg_m_s)
                .chain(m.angular_momentum_kg_m2_s)
                .all(f64::is_finite)
        };
        if !finite(&source)
            || !finite(&target)
            || !source_h.is_finite()
            || !target_h.is_finite()
            || (target.mass_kg - source.mass_kg).abs() > limits.mass_tolerance_kg
            || (liquid
                .flow()
                .ledger_mass_kg()
                .map_err(FreeSurfaceEnthalpyError::from)?
                - source.mass_kg)
                .abs()
                > limits.mass_tolerance_kg
            || (target_h - source_h).abs() > limits.enthalpy_tolerance_j
            || distance(source.momentum_kg_m_s, target.momentum_kg_m_s)
                > limits.momentum_tolerance_kg_m_s
            || distance(
                source.angular_momentum_kg_m2_s,
                target.angular_momentum_kg_m2_s,
            ) > limits.angular_momentum_tolerance_kg_m2_s
            || deficit < -limits.kinetic_roundoff_tolerance_j
            || deficit > limits.maximum_kinetic_projection_deficit_j
        {
            return Err(refused(
                "native liquid release exceeds its mass, heat, momentum or projection allowance",
            ));
        }
        checkpoint(cx, self.gate)?;
        Ok((
            liquid,
            LiquidReleaseReport {
                time_origin_s: thermal.time_s(),
                source,
                target,
                source_enthalpy_j: source_h,
                target_enthalpy_j: target_h,
                kinetic_projection_deficit_j: deficit,
            },
        ))
    }
}

fn checkpoint(cx: &Cx<'_>, gate: &CancelGate) -> Result<(), LiquidReleaseError> {
    if cx.checkpoint().is_err() || gate.is_requested() {
        Err(FreeSurfaceEnthalpyError::Flow(
            crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into(),
        )
        .into())
    } else {
        Ok(())
    }
}
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}
fn accumulate_motion(m: &mut LiquidTransferMoments, mass: f64, p: [f64; 3], v: [f64; 3]) {
    let momentum = v.map(|u| mass * u);
    for a in 0..3 {
        m.momentum_kg_m_s[a] += momentum[a];
    }
    for (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] {
        m.angular_momentum_kg_m2_s[a] += p[b] * momentum[c] - p[c] * momentum[b];
    }
    m.kinetic_energy_j += 0.5 * mass * v.iter().map(|u| u * u).sum::<f64>();
}
