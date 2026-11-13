//! Two finite aqueous inventories coupled to actual compatible swelling metrics.
//!
//! The original nested bath roots propose activities. Each activity pair enters
//! the original nonlinear diffusion/solid owner, which uses the candidate's
//! axial stretches and current surface/contact area. This closes geometric
//! feedback without creating a second diffusion, mechanics or phase model.

use super::*;
use crate::moisture::{
    layered::{MoistureActivityTransfer, atmosphere::water_change},
    reservoir::water_sum,
    solution::{
        FiniteMoistureSolution, MoistureSolutionPort, PairedSolutionMoistureCandidate,
        PairedSolutionMoistureSolve, PairedSolutionMoistureStep, SolutionMoistureError,
        advance_baths,
    },
    surface::admit_environment,
};

/// Both original bath roots and the shared native mechanical-work allowance.
#[derive(Clone, Copy, Debug)]
pub struct MoistureSpatialPairedSolutionControls {
    /// Original source-supported volume brackets, water tolerances and root
    /// limits. `maximum_evaluations` bounds all nested transport proposals.
    pub solutions: PairedSolutionMoistureSolve,
    /// Permitted departure from the separable axial prism map [m].
    pub maximum_geometry_residual_m: f64,
    /// All native solid solves across every bath and nonlinear-geometry trial,
    /// including rejected trials, in 1..=100_000. Reserved before each solve.
    pub maximum_native_solid_solves: usize,
}

/// Actual accepted metric and cumulative physical work of the coupled interval.
#[derive(Clone, Debug)]
pub struct MoistureSpatialPairedSolutionReport {
    /// Native source, transport, actual geometry and accepted-trial work.
    pub spatial: MoistureSpatialDiffusionReport,
    /// All native solves, including discarded bath and geometry proposals.
    /// The enclosing paired step counts all transport evaluations separately.
    pub native_solid_solves: usize,
}

/// Original finite solutions and layered water with swelling-dependent flux.
pub type MoistureSpatialPairedSolutionCandidate =
    PairedSolutionMoistureCandidate<MoistureLayeredSlab, MoistureSpatialPairedSolutionReport>;

impl MoistureLayeredSolid {
    /// Couple two independently finite baths to the actual swelling geometry.
    /// Each exterior face names First, Second, or None. A physically sealed
    /// transfer law remains sealed even if a bath is assigned; internal layer
    /// contacts never become reservoir boundaries. Each layer retains its own
    /// sourced sorption, diffusivity, dry inventory, elasticity and eigenstrain.
    ///
    /// For every trial, original native mechanics resolves compatible shape
    /// and the water residual uses that same candidate's inverse axial stretch
    /// squared and current surface/contact area. Only assigned signed face
    /// transfers debit each original bath. Its solute/phase owner is restarted
    /// from the original inventory for every root trial; no equilibrium reset
    /// replaces finite-rate dissolution or crystallization.
    ///
    /// This is first-order isothermal exchange through the equal-pressure
    /// vapor-equilibrium bridge and the existing separable-prism reduction.
    /// The scalar D is explicitly a CURRENT spatial diffusivity. Small-strain
    /// mechanics, source support and geometry admissibility remain hard gates.
    /// Shear, bending and taper outside the requested tolerance refuse. There
    /// is no immersed-liquid permeability, pressure-driven osmosis, inferred
    /// osmotic traction, stress-assisted chemical potential or thermal feedback.
    /// Specimen swelling does not displace either prescribed bath vessel.
    /// Chemical reservoir work is outside the reported mechanical energy.
    ///
    /// Both baths, the native moisture field and its actual solid state return
    /// together only after independent bath, field/face and whole-water audits.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_spatial_diffusion_paired_solutions(
        &self,
        solutions: [&FiniteMoistureSolution; 2],
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        faces: [Option<MoistureSolutionPort>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: MoistureSpatialPairedSolutionControls,
        cx: &Cx<'_>,
    ) -> Result<(MoistureSpatialPairedSolutionCandidate, MoistureSolidState), SolutionMoistureError>
    {
        self.advance_spatial_diffusion_paired_solutions_budgeted(
            solutions,
            moisture,
            dt_s,
            faces,
            transport,
            nonlinear,
            controls,
            &mut 0,
            || Ok(()),
            cx,
        )
    }

    /// The duration owner retains both counters across failed/rejected probes.
    /// Updating the mechanical counter before propagating a bath-root failure
    /// prevents a shortened retry from obtaining a fresh native-work allowance.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::moisture) fn advance_spatial_diffusion_paired_solutions_budgeted(
        &self,
        solutions: [&FiniteMoistureSolution; 2],
        moisture: &MoistureLayeredSlab,
        dt_s: f64,
        faces: [Option<MoistureSolutionPort>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        controls: MoistureSpatialPairedSolutionControls,
        native_solid_solves: &mut usize,
        charge_field: impl FnMut() -> Result<(), SolutionMoistureError>,
        cx: &Cx<'_>,
    ) -> Result<(MoistureSpatialPairedSolutionCandidate, MoistureSolidState), SolutionMoistureError>
    {
        poll(cx)?;
        if !positive(controls.maximum_geometry_residual_m)
            || !(1..=100_000).contains(&controls.maximum_native_solid_solves)
        {
            return Err(SolutionMoistureError::Invalid(
                "positive spatial geometry tolerance and bounded cumulative solid work required",
            ));
        }
        let assignments = std::array::from_fn(|side| {
            if moisture.exterior_transfer()[side] == MoistureActivityTransfer::Sealed {
                None
            } else {
                faces[side]
            }
        });
        let sealed = [
            !assignments.contains(&Some(MoistureSolutionPort::First)),
            !assignments.contains(&Some(MoistureSolutionPort::Second)),
        ];
        if assignments.iter().any(Option::is_some) {
            for layer in moisture.layers() {
                poll(cx)?;
                // The bridge supplies temperature and water activity only.
                // Additional pressure/composition coordinates need their owner.
                admit_environment(layer.sorption().source_point())?;
            }
        }
        let initial_work = *native_solid_solves;
        let mut work = SpatialSolveBudget {
            used: initial_work,
            maximum: controls.maximum_native_solid_solves,
        };
        let outcome = advance_baths(
            solutions,
            moisture.layers()[0].sorption().temperature_kelvin(),
            dt_s,
            sealed,
            controls.solutions,
            cx,
            |activities| {
                let boundary = assignments.map(|assignment| match assignment {
                    Some(port) => {
                        MoistureLayeredFace::PrescribedRelativeHumidity(activities[port.index()])
                    }
                    None => MoistureLayeredFace::Sealed,
                });
                let candidate = self.advance_spatial_diffusion_budgeted(
                    moisture,
                    dt_s,
                    &boundary,
                    transport,
                    nonlinear,
                    controls.maximum_geometry_residual_m,
                    &mut work,
                    cx,
                )?;
                let face_water = candidate.transport.transport.boundary_water_kg;
                let received =
                    std::array::from_fn(|port| {
                        water_sum(face_water.iter().zip(assignments).filter_map(
                            |(&water, face)| face.filter(|p| p.index() == port).map(|_| water),
                        ))
                    });
                Ok((candidate, received))
            },
            charge_field,
        );
        *native_solid_solves = work.used;
        let (baths, accepted, exchange, evaluations) = outcome?;
        let change = water_change(moisture, &accepted.moisture);
        let stored = |field: &MoistureLayeredSlab| {
            water_sum(field.layers().iter().flat_map(|layer| {
                let slab = layer.slab();
                slab.moisture()
                    .iter()
                    .enumerate()
                    .map(move |(i, m)| m * slab.field().dry_mass_at(i))
            }))
        };
        let solid_water_kg = [stored(moisture), stored(&accepted.moisture)];
        let field_balance = water_sum([
            change,
            -exchange[0].water_to_solid_kg,
            -exchange[1].water_to_solid_kg,
        ]);
        let bath_changes = [
            baths[0].solvent_mass_kg(),
            -solutions[0].solvent_mass_kg(),
            baths[1].solvent_mass_kg(),
            -solutions[1].solvent_mass_kg(),
        ];
        let balance = water_sum(bath_changes.into_iter().chain([change]));
        let direct = water_sum(
            bath_changes.into_iter().chain(
                moisture
                    .layers()
                    .iter()
                    .zip(accepted.moisture.layers())
                    .flat_map(|(a, b)| {
                        let old = a.slab();
                        let next = b.slab();
                        old.moisture()
                            .iter()
                            .zip(next.moisture())
                            .enumerate()
                            .flat_map(move |(i, (&a, &b))| {
                                [
                                    b * next.field().dry_mass_at(i),
                                    -a * old.field().dry_mass_at(i),
                                ]
                            })
                    }),
            ),
        );
        for residual_kg in [field_balance, balance, direct] {
            if !residual_kg.is_finite() || residual_kg.abs() > controls.solutions.mass_tolerance_kg
            {
                return Err(SolutionMoistureError::MassBalance {
                    residual_kg,
                    tolerance_kg: controls.solutions.mass_tolerance_kg,
                });
            }
        }
        if solid_water_kg.iter().any(|mass| !mass.is_finite()) {
            return Err(SolutionMoistureError::Invalid(
                "spatial specimen water inventory is unrepresentable",
            ));
        }
        let MoistureLayeredSolidStep {
            moisture,
            transport: spatial,
            solid,
        } = accepted;
        let report = MoistureSpatialPairedSolutionReport {
            spatial,
            native_solid_solves: work.used - initial_work,
        };
        let step = PairedSolutionMoistureStep {
            transport: report,
            face_solutions: assignments,
            solutions: exchange,
            solid_water_kg,
            water_change_kg: change,
            field_water_balance_residual_kg: field_balance,
            water_balance_residual_kg: balance,
            direct_water_balance_residual_kg: direct,
            field_evaluations: evaluations,
        };
        poll(cx)?;
        Ok((
            PairedSolutionMoistureCandidate::from_accepted(baths, moisture, step),
            solid,
        ))
    }
}
