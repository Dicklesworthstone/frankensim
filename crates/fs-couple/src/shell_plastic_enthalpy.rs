//! Reciprocal thermal plasticity and creep on the original layered CST/DKT shell.
//!
//! One original enthalpy cell owns each physical layer of each reference facet.
//! Temperature is uniform within that cell; all three surface points and every
//! thickness point retain independent native inelastic histories. The original
//! source families supply thermal expansion, softening, hardening, viscosity
//! and cyclic recovery, or standalone/joint primary and secondary creep. Their
//! physical and reversible heat enters the existing
//! thermal graph once. Numerical material loss and the retained reference
//! drilling gauge never become heat.
//!
//! Each physical interval has one endpoint equilibrium solve. Coupled trials
//! restart from the same accepted mechanical, thermal and observer histories;
//! only a balanced candidate is returned. This is fixed-reference, fully solid,
//! small-strain/rotation plane stress, not finite forming, transverse shear,
//! delamination or evolving thermal geometry. Refine physical intervals, native
//! thermal substeps, facet cells and thickness integration independently.

/// Sustained shell contact, bulk plastic heat and finite two-body heat transfer.
pub mod contact;
mod feedback;

use std::collections::BTreeSet;
use std::convert::Infallible;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::thermal::source::{
    ThermalSectionSourceError, ThermalSectionSourceFamily,
};
use fs_material::plastic::section::thermal::{ThermoplasticSectionError, ThermoplasticSectionStep};
use fs_material::plastic::section::{PlaneStressLayer, SectionBudget, ThicknessQuadrature};
use fs_material::plastic::thermal::ThermoplasticSource;
use fs_plate::shell::ShellMesh;
use fs_solid::plastic_shell::drilling::{PlasticShellReferenceDrilling, ShellDrillingEnergy};
use fs_solid::plastic_shell::dynamics::{
    PlasticShellDensity, PlasticShellReferenceInertia, ShellInertialEnergy,
};
use fs_solid::plastic_shell::{
    PlasticShellBudget, PlasticShellSections, SmallStrainShellProblem, SmallStrainShellSolution,
};
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionSimulation, EnergyAudit, ThermalLaw,
};

pub use crate::solid_enthalpy::plastic::ThermoplasticConfig;
pub use feedback::{
    ShellLayerHeat, ShellThermoplasticCandidate, ShellThermoplasticMotion, ShellThermoplasticReport,
};

type History<'a, S> = <<S as ThermalSectionSourceFamily<'a>>::Material as SmallStrainLaw>::State;

/// Original thermal identity and source-bound material layer, in increasing z.
#[derive(Clone, Debug)]
pub struct ShellThermoplasticLayer<S> {
    pub cell: CellId,
    pub layer: PlaneStressLayer<S>,
}

/// Immutable geometric reduction and explicit allocation limits.
#[derive(Clone, Copy, Debug)]
pub struct ShellThermoplasticDiscretization {
    pub coordinate_length_m: f64,
    pub rotation_limit_rad: f64,
    /// Fixed at binding; thermal softening cannot change this numerical gauge.
    pub drilling_alpha: f64,
    pub quadrature: ThicknessQuadrature,
    pub reduction: ReductionConfig,
    pub section_budget: SectionBudget,
    pub shell_budget: PlasticShellBudget,
}

/// Source, native section/shell, transport or coupled-balance refusal.
#[derive(Debug)]
pub enum ShellThermoplasticError<E = Infallible> {
    Invalid(&'static str),
    Material(ThermalSectionSourceError),
    Section(ThermoplasticSectionError),
    Solid(SmallStrainTetError),
    Contact(fs_solid::plastic_shell::contact::ShellContactError),
    Thermal(DepositionError),
    Observer(E),
    Cancelled,
    NotConverged {
        iterations: usize,
        temperature_defect_k: f64,
        heat_defect_j: f64,
        combined_residual_j: f64,
    },
}
impl<E: std::fmt::Debug> std::fmt::Display for ShellThermoplasticError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "shell thermoplastic coupling refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ShellThermoplasticError<E> {}
impl<E> From<DepositionError> for ShellThermoplasticError<E> {
    fn from(error: DepositionError) -> Self {
        if matches!(error, DepositionError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl<E> From<ThermoplasticSectionError> for ShellThermoplasticError<E> {
    fn from(error: ThermoplasticSectionError) -> Self {
        Self::Section(error)
    }
}

#[derive(Clone, Debug)]
struct Facet<S> {
    layers: Vec<PlaneStressLayer<S>>,
    reference_h_j_kg: Vec<f64>,
    area_m2: f64,
}

/// Accepted shell and immutable source/reference bindings. Native thermal
/// inventory is checked against `accepted`; it is not independently evolved.
#[derive(Clone, Debug)]
pub struct ThermoplasticShell<'a, S: ThermalSectionSourceFamily<'a> = ThermoplasticSource<'a>> {
    mesh: &'a ShellMesh,
    facets: Arc<[Facet<S>]>,
    accepted: Vec<Vec<DepositedCell>>,
    reference_drilling: Arc<PlasticShellReferenceDrilling<'a>>,
    solution: SmallStrainShellSolution<History<'a, S>>,
    discretization: ShellThermoplasticDiscretization,
    area_m2: f64,
    volume_m3: f64,
    time_s: f64,
}

impl<'a, S: ThermalSectionSourceFamily<'a>> ThermoplasticShell<'a, S> {
    /// Bind virgin, stress-free layers to their original active caloric objects.
    /// Sources may select J2, Voce, Perzyna, complete native Chaboche families,
    /// or `CaloricCreepSource` with original Norton/joint/expansion histories.
    /// Layer cells must be unique; mass is sourced density times reference
    /// facet area and physical layer thickness. No density or chart is invented.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        mesh: &'a ShellMesh,
        layers: Vec<Vec<ShellThermoplasticLayer<S>>>,
        discretization: ShellThermoplasticDiscretization,
        cx: &Cx<'_>,
    ) -> Result<Self, ShellThermoplasticError> {
        poll(cx)?;
        if layers.len() != mesh.tris.len()
            || layers.is_empty()
            || mesh.nodes.len() > discretization.shell_budget.maximum_nodes
            || layers.len() > discretization.shell_budget.maximum_elements
        {
            return Err(ShellThermoplasticError::Invalid(
                "bounded native mesh and one layer stack per facet required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut facets = Vec::with_capacity(layers.len());
        let mut accepted = Vec::with_capacity(layers.len());
        let mut area_m2 = 0.0;
        let mut volume_m3 = 0.0;
        let mut material_points = 0usize;
        for (e, entries) in layers.into_iter().enumerate() {
            poll(cx)?;
            if entries.is_empty() || entries.len() > discretization.section_budget.maximum_layers {
                return Err(ShellThermoplasticError::Invalid(
                    "bounded nonempty physical layer stack required",
                ));
            }
            let geometry = mesh
                .facet(e)
                .map_err(|_| ShellThermoplasticError::Invalid("invalid native shell facet"))?;
            let mut layer_sources = Vec::with_capacity(entries.len());
            let mut cells = Vec::with_capacity(entries.len());
            for entry in entries {
                poll(cx)?;
                let cell =
                    thermal
                        .cells()
                        .get(&entry.cell)
                        .ok_or(ShellThermoplasticError::Invalid(
                            "inactive shell layer cell",
                        ))?;
                let source = &entry.layer.law;
                let volume = geometry.area_m2 * (entry.layer.top_z_m - entry.layer.bottom_z_m);
                if !used.insert(entry.cell)
                    || !thermal.uses_material_law(
                        cell.material_id,
                        source.caloric_source().curve() as &dyn ThermalLaw,
                    )
                    || !same_measure(cell.geometry.volume_m3, volume)
                    || !same_measure(cell.mass_kg, source.reference_density_kg_m3() * volume)
                {
                    return Err(ShellThermoplasticError::Invalid(
                        "unique original layer cells, caloric objects and sourced reference inventory required",
                    ));
                }
                // An already heated source with eigenstrain needs an explicit
                // initial equilibrium, not a fabricated stress-free checkpoint.
                let step = entry
                    .layer
                    .law
                    .resolve_section_step(cell.specific_h_j_kg, cell.specific_h_j_kg, 1.0)
                    .map_err(|source| ThermoplasticSectionError::Source {
                        layer: cells.len(),
                        point: None,
                        source,
                    })?;
                let virgin = step.law.initial_state();
                if step
                    .law
                    .stress(&[0.0; 6], &virgin)
                    .iter()
                    .any(|&x| x != 0.0)
                    || step.law.update_state(&[0.0; 6], &virgin) != virgin
                {
                    return Err(ShellThermoplasticError::Invalid(
                        "virgin sources must be stress free at the actual initial enthalpy",
                    ));
                }
                material_points = entry
                    .layer
                    .subdivisions
                    .checked_mul(discretization.quadrature.points_per_subdivision())
                    .and_then(|n| n.checked_mul(3))
                    .and_then(|n| material_points.checked_add(n))
                    .ok_or(ShellThermoplasticError::Invalid(
                        "shell history count overflow",
                    ))?;
                if material_points > discretization.shell_budget.maximum_material_points {
                    return Err(ShellThermoplasticError::Invalid(
                        "shell material history budget exceeded",
                    ));
                }
                volume_m3 = finite(volume_m3 + volume)?;
                layer_sources.push(entry.layer);
                cells.push(cell.clone());
            }
            area_m2 = finite(area_m2 + geometry.area_m2)?;
            facets.push(Facet {
                layers: layer_sources,
                reference_h_j_kg: cells.iter().map(|cell| cell.specific_h_j_kg).collect(),
                area_m2: geometry.area_m2,
            });
            accepted.push(cells);
        }
        let enthalpies: Vec<Vec<_>> = accepted
            .iter()
            .map(|cells| cells.iter().map(|cell| cell.specific_h_j_kg).collect())
            .collect();
        let sections: Vec<_> = facets
            .iter()
            .zip(&enthalpies)
            .map(|(facet, h)| {
                ThermoplasticSectionStep::new(
                    &facet.layers,
                    h,
                    h,
                    1.0,
                    discretization.quadrature,
                    discretization.reduction,
                    discretization.section_budget,
                )
            })
            .collect::<Result<_, _>>()?;
        let references: Vec<_> = sections
            .iter()
            .map(ThermoplasticSectionStep::section)
            .collect();
        // Reference preparation performs no equilibrium solve. Constrain the
        // virgin coordinates so it does not consume the later free-DOF budget.
        let reference_supports: Vec<_> = (0..6 * mesh.nodes.len()).map(|d| (d, 0.0)).collect();
        let problem = SmallStrainShellProblem {
            mesh,
            sections: PlasticShellSections::PerElementBorrowed(&references),
            coordinate_length_m: discretization.coordinate_length_m,
            rotation_limit_rad: discretization.rotation_limit_rad,
            drilling_alpha: discretization.drilling_alpha,
            prescribed: &reference_supports,
            nodal_loads: &[],
            load_path: &[1.0],
            budget: discretization.shell_budget,
            settings: IncrementSettings::default(),
        };
        let reference_drilling = Arc::new(
            PlasticShellReferenceDrilling::prepare(mesh, &problem, cx)
                .map_err(ShellThermoplasticError::Solid)?,
        );
        let solution = SmallStrainShellSolution {
            coordinates: vec![[0.0; 6]; mesh.nodes.len()],
            section_states: sections
                .iter()
                .flat_map(|section| (0..3).map(move |_| section.section().initial_state()))
                .collect(),
            section_resultants: vec![[0.0; 6]; 3 * facets.len()],
            increments: Vec::new(),
        };
        poll(cx)?;
        Ok(Self {
            mesh,
            facets: facets.into(),
            accepted,
            reference_drilling,
            solution,
            discretization,
            area_m2,
            volume_m3,
            time_s: thermal.time_s(),
        })
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn solution(&self) -> &SmallStrainShellSolution<History<'a, S>> {
        &self.solution
    }

    /// Prepare the original consistent physical/director mass once for repeated
    /// dynamic intervals. Density comes from the admitted source inventory;
    /// later heating cannot create mass. Static trajectories need not allocate
    /// this operator. It borrows the original mesh, not temporary thermal laws.
    pub fn prepare_reference_inertia(
        &self,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellReferenceInertia<'a>, ShellThermoplasticError> {
        poll(cx)?;
        let d = self.discretization;
        let mut sections = Vec::with_capacity(self.facets.len());
        let mut density = Vec::with_capacity(self.facets.len());
        for facet in self.facets.iter() {
            poll(cx)?;
            sections.push(ThermoplasticSectionStep::new(
                &facet.layers,
                &facet.reference_h_j_kg,
                &facet.reference_h_j_kg,
                1.0,
                d.quadrature,
                d.reduction,
                d.section_budget,
            )?);
            density.push(
                facet
                    .layers
                    .iter()
                    .map(|layer| layer.law.reference_density_kg_m3())
                    .collect(),
            );
        }
        let references: Vec<_> = sections
            .iter()
            .map(ThermoplasticSectionStep::section)
            .collect();
        let reference_supports: Vec<_> = (0..6 * self.mesh.nodes.len()).map(|d| (d, 0.0)).collect();
        let problem = SmallStrainShellProblem {
            mesh: self.mesh,
            sections: PlasticShellSections::PerElementBorrowed(&references),
            coordinate_length_m: d.coordinate_length_m,
            rotation_limit_rad: d.rotation_limit_rad,
            drilling_alpha: d.drilling_alpha,
            prescribed: &reference_supports,
            nodal_loads: &[],
            load_path: &[1.0],
            budget: d.shell_budget,
            settings: IncrementSettings::default(),
        };
        PlasticShellReferenceInertia::prepare(
            self.mesh,
            &problem,
            PlasticShellDensity::PerElement(&density),
            cx,
        )
        .map_err(ShellThermoplasticError::Solid)
    }

    /// Resolve transport from these same source identities. Conductivity is
    /// explicitly selected per physical layer in facet order; the existing
    /// material-network observer owns conduction and its thermal substeps.
    pub fn thermal_material_cells<'s>(
        &'s self,
        conductivity: &'s [Option<(
            fs_matdb::PropertyKey,
            fs_material::state_point::MaterialPropertySelection,
        )>],
    ) -> Result<Vec<crate::circuit_enthalpy::material::MaterialCell<'s>>, ShellThermoplasticError>
    {
        if conductivity.len() != self.accepted.iter().map(Vec::len).sum::<usize>() {
            return Err(ShellThermoplasticError::Invalid(
                "one transport selection per shell layer cell required",
            ));
        }
        self.facets
            .iter()
            .zip(&self.accepted)
            .flat_map(|(facet, cells)| facet.layers.iter().zip(cells))
            .zip(conductivity)
            .map(|((layer, cell), selection)| {
                let source = &layer.law;
                let conductivity = selection
                    .as_ref()
                    .map(|(key, policy)| {
                        source
                            .pin_thermal_transport(cell.specific_h_j_kg, key, policy.clone())
                            .map(|pinned| (key, pinned))
                            .map_err(ShellThermoplasticError::Material)
                    })
                    .transpose()?;
                Ok(crate::circuit_enthalpy::material::MaterialCell {
                    cell: cell.id,
                    material_id: cell.material_id,
                    curve: source.caloric_source().curve(),
                    point: source.thermal_material_point(),
                    conductivity,
                })
            })
            .collect()
    }

    fn enthalpies<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        require_accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<Vec<Vec<f64>>, ShellThermoplasticError<E>> {
        let mut result = Vec::with_capacity(self.facets.len());
        for (facet, cells) in self.facets.iter().zip(&self.accepted) {
            let mut h = Vec::with_capacity(cells.len());
            for (layer, old) in facet.layers.iter().zip(cells) {
                poll(cx)?;
                let cell = thermal
                    .cells()
                    .get(&old.id)
                    .ok_or(ShellThermoplasticError::Invalid(
                        "shell layer cell became inactive",
                    ))?;
                let source = &layer.law;
                if cell.id != old.id
                    || cell.segment_id != old.segment_id
                    || cell.body_id != old.body_id
                    || cell.material_id != old.material_id
                    || cell.geometry != old.geometry
                    || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
                    || require_accepted
                        && cell.specific_h_j_kg.to_bits() != old.specific_h_j_kg.to_bits()
                    || !thermal.uses_material_law(
                        cell.material_id,
                        source.caloric_source().curve() as &dyn ThermalLaw,
                    )
                {
                    return Err(ShellThermoplasticError::Invalid(
                        "shell thermal identity, reference inventory or accepted enthalpy changed",
                    ));
                }
                source
                    .admit_section_enthalpy(cell.specific_h_j_kg)
                    .map_err(ShellThermoplasticError::Material)?;
                h.push(cell.specific_h_j_kg);
            }
            result.push(h);
        }
        Ok(result)
    }
}

/// Whole-shell physical and numerical accounts in joules, integrated over the
/// original facet area and exact retained thickness quadrature.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ShellThermoplasticEnergy {
    pub physical_heat_j: f64,
    /// Physical creep work, included exactly once in `physical_heat_j`.
    pub creep_heat_j: f64,
    pub reversible_heat_j: f64,
    /// Free-energy change including the separately reported coefficient change.
    pub stored_change_j: f64,
    pub coefficient_storage_change_j: f64,
    pub thermoelastic_storage_change_j: f64,
    pub material_numerical_loss_j: f64,
    pub material_work_j: f64,
    pub generalized_work_j: f64,
    pub transverse_residual_work_j: f64,
    pub drilling: ShellDrillingEnergy,
    pub dead_load_work_j: f64,
    pub support_work_j: f64,
    /// Work of actual solved normal and tangential forces on the shell.
    /// A coupled interface replaces this port by obstacle work and the distinct
    /// physical friction/numerical approach accounts; it is not all heat.
    pub contact_work_j: f64,
    /// Signed endpoint inertia work; physical and gauge kinetics are separately
    /// exposed in the dynamic candidate. No inertial numerical loss is heat.
    pub inertial_work_j: f64,
    pub equilibrium_residual_j: f64,
}
impl ShellThermoplasticEnergy {
    pub fn external_work_j(&self) -> f64 {
        self.dead_load_work_j + self.support_work_j + self.contact_work_j
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), ShellThermoplasticError<E>> {
    cx.checkpoint()
        .map_err(|_| ShellThermoplasticError::Cancelled)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite<E>(value: f64) -> Result<f64, ShellThermoplasticError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ShellThermoplasticError::Invalid(
            "unrepresentable shell heat or work arithmetic",
        ))
    }
}
fn same_measure(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 512.0 * f64::EPSILON * a.abs().max(b.abs())
}
