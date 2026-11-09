//! Source-bound saturating inductors on the original Kirchhoff descriptor.
//!
//! A uniform core keeps every original H(B) breakpoint. With core flux Phi=A B,
//! N I = l H(B) + R_gap Phi and lambda = N Phi + L_leak I. Both mappings are
//! piecewise affine, so their composition gives the exact native I(lambda)
//! curve without sampling, fitting, extrapolating or solving a second ODE.
//! Its integral is core energy + gap energy + leakage energy:
//! A*l*integral H dB + R_gap*Phi^2/2 + L_leak*I^2/2.
//!
//! This is an explicitly declared uniform-flux magnetic-path reduction.
//! Geometry, source context, gap reluctance and leakage stay fixed. Fringing,
//! distributed leakage, hysteresis, eddy loss and thermal magnetic entropy need
//! their own physical owners. Graph resistance and other storage remain native.
//! Energy/flux conventions follow MIT 6.685 Class Notes 1:
//! <https://web.mit.edu/6.685/www/chapter1.pdf>.

/// Source-bound multiwinding saturation on one common magnetic path.
pub mod shared;
pub use shared::{
    MaterialSharedCoreInductor, ReversibleCoreGeometry, SharedCoreOperatingPoint, SharedCoreWinding,
};

use fs_material::field_storage::{
    nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
    permanent::FLUX_DENSITY_DIMS,
};
use fs_phs::circuit::{
    Branch, CircuitDae, CircuitGraph, MutualInductance,
    capacitance::{CapacitanceBlock, CapacitorCurve},
    inductance::{
        InductorCurve, InductorFluxCurve, MAX_INDUCTOR_KNOTS, assemble_circuit_with_storage_blocks,
    },
};
use std::sync::Arc;

use super::{
    MAX_MATERIAL_CIRCUIT_BRANCHES, MAX_MATERIAL_CIRCUIT_NODES, MaterialCircuitError, admit_graph,
    extended, poll, positive,
};

fn invalid(message: &'static str) -> MaterialCircuitError {
    MaterialCircuitError::Invalid(message)
}

fn finite(value: f64) -> Result<f64, MaterialCircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "reversible inductor arithmetic is not representable",
        ))
    }
}

/// Fixed, explicitly authored one-dimensional magnetic-path geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReversibleInductorGeometry {
    /// Positive effective winding turns; current orientation follows its branch.
    pub turns: f64,
    /// Uniform core cross-section [m^2].
    pub core_area_m2: f64,
    /// Uniform material path length [m].
    pub core_length_m: f64,
    /// Additional linear series reluctance [1/H], e.g. gap/(mu_gap*area).
    /// Zero omits that region. No fringing or air permeability is inferred.
    pub gap_reluctance_per_h: f64,
    /// Explicit linear series leakage inductance [H], outside the core/gap path.
    pub leakage_inductance_h: f64,
}

/// Immutable source, geometry and native branch law, sharing prepared tables.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialReversibleInductor {
    branch: usize,
    material: ReversibleMagneticMaterial,
    geometry: ReversibleInductorGeometry,
    curve: InductorFluxCurve,
    flux_density_at_knots_t: Arc<[f64]>,
}

/// Physical core/gap/leakage observations at one admitted flux linkage.
/// At a source kink, differential inductance uses the native curve's right-linkage
/// derivative (left derivative only at its final endpoint).
/// Energy components partition the native stored energy; they are not extra
/// stores to add to the enclosing circuit a second time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReversibleInductorOperatingPoint {
    pub branch: usize,
    pub flux_linkage_wb: f64,
    pub current_a: f64,
    pub core_flux_density_t: f64,
    pub core_field_strength_a_m: f64,
    pub differential_inductance_h: f64,
    pub core_energy_j: f64,
    pub gap_energy_j: f64,
    pub leakage_energy_j: f64,
    pub stored_energy_j: f64,
}

impl MaterialReversibleInductor {
    /// Map the whole admitted source path to both polarities of the isotropic
    /// zero-remanence law. Every physical breakpoint is retained exactly.
    ///
    /// A field-only source potential is required. Actual dependence on another
    /// response axis or an explicit thermal scale refuses: freezing a thermal
    /// free-energy snapshot cannot account for its magnetic entropy exchange.
    /// Fixed source applicability bounds and query context remain attached.
    pub fn new(
        branch: usize,
        material: ReversibleMagneticMaterial,
        geometry: ReversibleInductorGeometry,
    ) -> Result<Self, MaterialCircuitError> {
        let g = geometry;
        if branch >= MAX_MATERIAL_CIRCUIT_BRANCHES
            || !positive(g.turns)
            || !positive(g.core_area_m2)
            || !positive(g.core_length_m)
            || !positive(g.core_area_m2 * g.core_length_m)
            || !positive(g.turns * g.core_area_m2)
            || !g.gap_reluctance_per_h.is_finite()
            || g.gap_reluctance_per_h < 0.0
            || !g.leakage_inductance_h.is_finite()
            || g.leakage_inductance_h < 0.0
        {
            return Err(invalid(
                "bounded inductor branch and finite positive core geometry required",
            ));
        }
        let property = material
            .material()
            .property(REVERSIBLE_FIELD_STRENGTH_PROPERTY)
            .ok_or_else(|| invalid("nonlinear inductor lost its original magnetic source"))?;
        let axes = property.source_axes();
        if axes.len() != 1
            || axes[0].0 != material.source_path().axis
            || axes[0].1 != FLUX_DENSITY_DIMS
            || material.field_scale_slope_per_k().is_some()
        {
            return Err(invalid(
                "closed nonlinear inductor requires a field-only source potential",
            ));
        }
        let path = &material.source_path().knots;
        let count = path
            .len()
            .checked_mul(2)
            .and_then(|n| n.checked_sub(1))
            .filter(|&n| n <= MAX_INDUCTOR_KNOTS)
            .ok_or_else(|| invalid("magnetic source path exceeds signed inductor knot capacity"))?;
        let mut knots = Vec::with_capacity(count);
        let mut fields = Vec::with_capacity(count);
        let transform = |b: f64| -> Result<(f64, f64), MaterialCircuitError> {
            let h = material.evaluate(b)?.field_strength_a_m;
            let phi = finite(g.core_area_m2 * b)?;
            let current =
                finite(finite(g.core_length_m * h + g.gap_reluctance_per_h * phi)? / g.turns)?;
            let linkage = finite(g.leakage_inductance_h.mul_add(current, g.turns * phi))?;
            Ok((linkage, current))
        };
        for &(b, _) in path.iter().skip(1).rev() {
            let (linkage, current) = transform(b)?;
            knots.push((-linkage, -current));
            fields.push(-b);
        }
        for &(b, _) in path {
            knots.push(transform(b)?);
            fields.push(b);
        }
        let result = Self {
            branch,
            material,
            geometry,
            curve: InductorFluxCurve::try_new(knots)?,
            flux_density_at_knots_t: fields.into(),
        };
        // Nonnegative energy components attain their largest value at the
        // endpoint. Refuse a reduction whose physical observations overflow.
        result.operating_point(result.curve.flux_range_wb().1)?;
        Ok(result)
    }

    pub const fn branch(&self) -> usize {
        self.branch
    }
    pub fn material(&self) -> &ReversibleMagneticMaterial {
        &self.material
    }
    pub const fn geometry(&self) -> ReversibleInductorGeometry {
        self.geometry
    }
    pub fn curve(&self) -> &InductorFluxCurve {
        &self.curve
    }
    /// Exact native branch binding, with immutable tables shared on clone.
    pub fn binding(&self) -> InductorCurve {
        InductorCurve {
            branch: self.branch,
            curve: self.curve.clone(),
        }
    }
    pub fn current_a(&self, flux_linkage_wb: f64) -> Result<f64, MaterialCircuitError> {
        Ok(self.curve.current_a(flux_linkage_wb)?)
    }
    pub fn flux_linkage_wb(&self, current_a: f64) -> Result<f64, MaterialCircuitError> {
        Ok(self.curve.flux_wb(current_a)?)
    }
    pub fn stored_energy_j(&self, flux_linkage_wb: f64) -> Result<f64, MaterialCircuitError> {
        Ok(self.curve.stored_energy_j(flux_linkage_wb)?)
    }

    /// Recover physical core fields using the original transformed segments.
    /// Interpolating B on lambda avoids cancellation of lambda-L_leak*I when
    /// declared leakage dominates the total linkage.
    pub fn operating_point(
        &self,
        flux_linkage_wb: f64,
    ) -> Result<ReversibleInductorOperatingPoint, MaterialCircuitError> {
        let current = self.curve.current_a(flux_linkage_wb)?;
        let knots = self.curve.knots();
        let segment = knots
            .partition_point(|&(lambda, _)| lambda <= flux_linkage_wb)
            .saturating_sub(1)
            .min(knots.len() - 2);
        let anchor = if flux_linkage_wb < 0.0 {
            segment + 1
        } else {
            segment
        };
        let b = if flux_linkage_wb == knots[segment].0 {
            self.flux_density_at_knots_t[segment]
        } else if flux_linkage_wb == knots[segment + 1].0 {
            self.flux_density_at_knots_t[segment + 1]
        } else {
            let fraction = finite(
                (flux_linkage_wb - knots[anchor].0) / (knots[segment + 1].0 - knots[segment].0),
            )?;
            finite(fraction.mul_add(
                self.flux_density_at_knots_t[segment + 1] - self.flux_density_at_knots_t[segment],
                self.flux_density_at_knots_t[anchor],
            ))?
        };
        let response = self.material.evaluate(b.abs())?;
        let g = self.geometry;
        let phi = finite(g.core_area_m2 * b)?;
        let core_energy_j =
            finite((g.core_area_m2 * g.core_length_m) * response.energy_density_j_m3)?;
        let gap_energy_j = finite((0.5 * g.gap_reluctance_per_h * phi) * phi)?;
        let leakage_energy_j = finite((0.5 * g.leakage_inductance_h * current) * current)?;
        Ok(ReversibleInductorOperatingPoint {
            branch: self.branch,
            flux_linkage_wb,
            current_a: current,
            core_flux_density_t: b,
            core_field_strength_a_m: b.signum() * response.field_strength_a_m,
            differential_inductance_h: self.curve.differential_inductance_h(flux_linkage_wb)?,
            core_energy_j,
            gap_energy_j,
            leakage_energy_j,
            stored_energy_j: self.curve.stored_energy_j(flux_linkage_wb)?,
        })
    }
}

/// Native circuit plus retained source bindings. This immutable assembly owns
/// no additional time history or integrator: use its original CircuitDae for
/// stepping, source work and resistor dissipation.
pub struct ReversibleInductorCircuit {
    graph: CircuitGraph,
    dae: CircuitDae,
    sources: Vec<MaterialReversibleInductor>,
    flux_indices: Vec<usize>,
    shared_cores: Vec<MaterialSharedCoreInductor>,
    shared_flux_indices: Vec<Vec<usize>>,
}

impl ReversibleInductorCircuit {
    /// Compose sourced saturating inductors with authored scalar storage,
    /// reciprocal capacitance blocks, and disjoint linear mutual inductance.
    /// Source-bound branch coefficients are replaced by their zero-field
    /// differential L before native topology admission; they are not fallbacks.
    ///
    /// Every original source remains inspectable. The native builder rejects
    /// duplicate branch ownership and any mutual link touching a scalar
    /// nonlinear inductor. Cancellation refuses before publishing the assembly.
    pub fn new(
        graph: CircuitGraph,
        sources: Vec<MaterialReversibleInductor>,
        capacitance_blocks: &[CapacitanceBlock],
        capacitor_curves: &[CapacitorCurve],
        authored_inductor_curves: &[InductorCurve],
        mutual: &[MutualInductance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        Self::with_shared_cores(
            graph,
            sources,
            vec![],
            capacitance_blocks,
            capacitor_curves,
            authored_inductor_curves,
            mutual,
            cancelled,
        )
    }

    /// Compose source-bound shared cores with the same scalar storage and
    /// circuit topology. Each shared core owns ALL self and mutual storage of
    /// its declared windings, including their explicit physical leakage.
    /// Branch ownership must be disjoint from other cores, scalar laws and
    /// authored linear mutual links. Original flux states and native solver
    /// remain authoritative; no private winding ports or outer solve are used.
    #[allow(clippy::too_many_arguments)]
    pub fn with_shared_cores(
        mut graph: CircuitGraph,
        mut sources: Vec<MaterialReversibleInductor>,
        mut shared_cores: Vec<MaterialSharedCoreInductor>,
        capacitance_blocks: &[CapacitanceBlock],
        capacitor_curves: &[CapacitorCurve],
        authored_inductor_curves: &[InductorCurve],
        mutual: &[MutualInductance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        poll(cancelled)?;
        if graph.node_count > MAX_MATERIAL_CIRCUIT_NODES
            || graph.branches.len() > MAX_MATERIAL_CIRCUIT_BRANCHES
            || sources
                .len()
                .checked_add(authored_inductor_curves.len())
                .and_then(|count| {
                    shared_cores
                        .iter()
                        .try_fold(count, |n, core| n.checked_add(core.windings().len()))
                })
                .is_none_or(|count| count > graph.branches.len())
            || capacitor_curves.len() > graph.branches.len()
            || capacitance_blocks
                .iter()
                .try_fold(0usize, |count, block| {
                    count.checked_add(block.branches().len())
                })
                .is_none_or(|count| count > graph.branches.len())
            || mutual.len()
                > graph
                    .branches
                    .len()
                    .saturating_mul(graph.branches.len().saturating_sub(1))
                    / 2
        {
            return Err(invalid(
                "source-bound inductor circuit exceeds native graph capacity",
            ));
        }
        sources.sort_by_key(MaterialReversibleInductor::branch);
        if sources
            .windows(2)
            .any(|pair| pair[0].branch == pair[1].branch)
        {
            return Err(invalid("duplicate source-bound inductor branch"));
        }
        let mut curves = Vec::with_capacity(sources.len() + authored_inductor_curves.len());
        for source in &sources {
            poll(cancelled)?;
            let branch = graph
                .branches
                .get_mut(source.branch)
                .ok_or_else(|| invalid("source-bound inductor branch is absent"))?;
            if !matches!(branch.2, Branch::Inductor { .. }) {
                return Err(invalid(
                    "magnetic source requires an original inductor branch",
                ));
            }
            branch.2 = Branch::Inductor {
                henries: source.curve.differential_inductance_h(0.0)?,
            };
            curves.push(source.binding());
        }
        curves.extend_from_slice(authored_inductor_curves);
        // Canonical core order does not change the signed winding order within
        // a core, which belongs to its authored magnetic geometry.
        shared_cores.sort_by_key(|core| core.windings().iter().map(|winding| winding.branch).min());
        let mut magnetic_blocks = Vec::with_capacity(shared_cores.len());
        for core in &shared_cores {
            poll(cancelled)?;
            let count = core.windings().len();
            let tangent = core.law().differential_inductance_h(&vec![0.0; count])?;
            for (i, winding) in core.windings().iter().enumerate() {
                let branch = graph
                    .branches
                    .get_mut(winding.branch)
                    .ok_or_else(|| invalid("shared-core winding branch is absent"))?;
                if !matches!(branch.2, Branch::Inductor { .. }) {
                    return Err(invalid("shared core requires original inductor branches"));
                }
                branch.2 = Branch::Inductor {
                    henries: tangent[i * count + i],
                };
            }
            magnetic_blocks.push(core.binding());
        }
        admit_graph(&graph)?;
        let dae = assemble_circuit_with_storage_blocks(
            &graph,
            capacitance_blocks,
            capacitor_curves,
            &curves,
            &magnetic_blocks,
            mutual,
        )?;
        let mut flux_indices = Vec::with_capacity(sources.len());
        for source in &sources {
            poll(cancelled)?;
            let slot = graph.branches[..source.branch]
                .iter()
                .filter(|branch| matches!(branch.2, Branch::Inductor { .. }))
                .count();
            flux_indices.push(dae.flux_index[slot]);
        }
        let mut shared_flux_indices = Vec::with_capacity(shared_cores.len());
        for core in &shared_cores {
            poll(cancelled)?;
            shared_flux_indices.push(
                core.windings()
                    .iter()
                    .map(|winding| {
                        let slot = graph.branches[..winding.branch]
                            .iter()
                            .filter(|branch| matches!(branch.2, Branch::Inductor { .. }))
                            .count();
                        dae.flux_index[slot]
                    })
                    .collect(),
            );
        }
        poll(cancelled)?;
        Ok(Self {
            graph,
            dae,
            sources,
            flux_indices,
            shared_cores,
            shared_flux_indices,
        })
    }

    pub fn graph(&self) -> &CircuitGraph {
        &self.graph
    }
    /// Original Kirchhoff equations, analytic constitutive derivatives and
    /// energy/work-audited time step. No separate magnetic solver is introduced.
    pub fn dae(&self) -> &CircuitDae {
        &self.dae
    }
    pub fn sources(&self) -> &[MaterialReversibleInductor] {
        &self.sources
    }

    /// Common material/geometry owners, ordered by their lowest branch index.
    pub fn shared_cores(&self) -> &[MaterialSharedCoreInductor] {
        &self.shared_cores
    }

    /// Initialize original branch-ordered fluxes and charges without changing
    /// those conserved coordinates to satisfy an incompatible source.
    pub fn consistent_initial_state(
        &self,
        fluxes_wb: &[f64],
        charges_c: &[f64],
        inputs: &[f64],
    ) -> Result<Vec<f64>, MaterialCircuitError> {
        let initial = extended(&self.dae, fluxes_wb, charges_c, inputs)?;
        Ok(self.dae.consistent_initial_state(&initial, inputs)?)
    }

    /// Source/core observations in stable original branch order.
    pub fn operating_points(
        &self,
        state: &[f64],
    ) -> Result<Vec<ReversibleInductorOperatingPoint>, MaterialCircuitError> {
        if state.len() != self.dae.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite native circuit state required"));
        }
        self.sources
            .iter()
            .zip(&self.flux_indices)
            .map(|(source, &index)| source.operating_point(state[index]))
            .collect()
    }

    /// One observation per common core; its core/gap energy is counted once,
    /// with separate leakage energy for every signed winding.
    pub fn shared_operating_points(
        &self,
        state: &[f64],
    ) -> Result<Vec<SharedCoreOperatingPoint>, MaterialCircuitError> {
        if state.len() != self.dae.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite native circuit state required"));
        }
        self.shared_cores
            .iter()
            .zip(&self.shared_flux_indices)
            .map(|(core, indices)| {
                let fluxes: Vec<_> = indices.iter().map(|&index| state[index]).collect();
                core.operating_point(&fluxes)
            })
            .collect()
    }
}
