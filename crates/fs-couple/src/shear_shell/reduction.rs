//! One elastic eigenbasis for shell force, volume and full material loss.
//! Discretization uses the original MITC3+ owner and fs-modal spectral slice.
//! No off-diagonal damping coefficient is discarded; pressure and cavity motion
//! observe the same physical displacement vectors as every reciprocal bridge.

use super::{invalid, owner, poll, ResolvedShearShell};
use crate::acoustic_realize::AcousticRealizeError;
use crate::plate_viscosity::ModalPlateViscosity;
use crate::thin_plate::CompactBody;
use fs_exec::CancelGate;
use fs_plate::shell::{modes_shell, mitc::{assemble_mitc3p, MitcBudget, NodalDirectors,
    ShellConstraint, ShellMass}};

/// Largest retained mode count one shell body or shell assembly reduction admits. It bounds the
/// stored mode vectors and the dense modal damping factor; the cold projection work still has to fit
/// the caller's `max_projection_terms`, which is the cost governor, and the mesh must actually hold
/// that many modes inside its eigenvalue window.
pub const MAX_SHELL_MODES: usize = 256;
// The shell connection projection inside fs-plate admits at most this many physical modes; a larger ceiling here
// would pass the check below and fail deep inside the eigen reduction.
const _: () = assert!(MAX_SHELL_MODES <= fs_plate::shell::mitc::connections::MAX_CONNECTION_MODES);

/// Whether `n` retained modes fit a declared `max_modes` that is itself within
/// `1..=`[`MAX_SHELL_MODES`]. Shared by the single-body and the assembly reductions.
#[must_use]
pub const fn shell_mode_budget_admits(n: usize, max_modes: usize) -> bool {
    n >= 1 && n <= max_modes && max_modes >= 1 && max_modes <= MAX_SHELL_MODES
}

/// Explicit cold projection resources and passive-channel coordinate scaling.
#[derive(Clone, Copy, Debug)]
pub struct ShellProjectionBudget {
    /// Complete requested mode count must fit (1..=[`MAX_SHELL_MODES`]).
    pub max_modes: usize,
    /// Condensation, all mode/force/volume projections and cubic factor work.
    pub max_projection_terms: usize,
    /// Positive dashpot resistance scale [N s/m]; C is unchanged by this scale.
    pub reference_resistance_n_s_m: f64,
}
/// Geometry-independent controls for the existing shell assembly and eigenowner.
#[derive(Clone, Debug)]
pub struct ShellReductionOptions {
    /// Actual homogeneous support constraints in physical or local axes.
    pub constraints: Vec<ShellConstraint>,
    /// Explicit translational and rotary mass discretization.
    pub mass: ShellMass,
    /// Positive drilling stabilization; contributes no fictitious material loss.
    pub drilling_alpha: f64,
    /// Angular-frequency-squared interval (low,high] [rad²/s²].
    pub eigenvalue_window: (f64,f64),
    /// Retain the lowest n modes inside the declared window, with no refit.
    pub n_modes: usize,
    /// Explicit assembly extent.
    pub assembly_budget: MitcBudget,
    /// Complete cold force/loss/volume projection allowance.
    pub projection_budget: ShellProjectionBudget,
}
/// Immutable common elastic basis and every work-conjugate projection.
#[derive(Clone, Debug)]
pub struct MaterialShellImage {
    source: ResolvedShearShell,
    directors: NodalDirectors,
    physical_modes: Vec<Vec<f64>>,
    force_ports: Vec<Vec<f64>>,
    volume_shapes: Vec<f64>,
    cavity_volume_shapes: Vec<f64>,
    intervals: Vec<(f64,f64)>,
    radiators: Vec<CompactBody>,
    damping: ModalPlateViscosity,
    interface_columns: Vec<Vec<f64>>,
    coating_relaxation: Option<super::coating::relaxation::ShellCoatingRelaxation>,
    substrate_relaxation: Option<super::constant_loss::relaxation::ShellSubstrateRelaxation>,
}
impl MaterialShellImage {
    /// Exact original material, section and midsurface sources.
    pub const fn source(&self) -> &ResolvedShearShell { &self.source }
    /// Exact original directors used by the MITC assembly and bubble condensation.
    pub const fn directors(&self) -> &NodalDirectors { &self.directors }
    /// Mode-major GLOBAL displacements/rotations, six entries per physical node.
    /// Vectors retain the original fs-modal mass normalization and mode signs.
    pub fn physical_modes(&self) -> &[Vec<f64>] { &self.physical_modes }
    /// Work-conjugate modal force/displacement rows in original caller order.
    pub fn force_ports(&self) -> &[Vec<f64>] { &self.force_ports }
    /// Exterior integral n dot phi dA, [m²/sqrt(kg)] in the same normalized basis.
    /// Only the explicitly selected exterior facets enter compact observation.
    /// Zero monopole volume does not establish silence under other radiation.
    pub fn volume_shapes_m2_per_sqrt_kg(&self) -> &[f64] { &self.volume_shapes }
    /// Independently selected outward-from-cavity integral n dot phi dA in the
    /// SAME basis [m²/sqrt(kg)]. Compression work uses this row once, regardless
    /// of exterior observation or the number of internal structural members.
    pub fn cavity_volume_shapes_m2_per_sqrt_kg(&self) -> &[f64] { &self.cavity_volume_shapes }
    /// Original eigenvalue residual intervals, not continuum error bounds.
    pub fn eigenvalue_intervals(&self) -> &[(f64,f64)] { &self.intervals }
    /// Complete symmetric modal C and passive connection realization.
    pub const fn damping(&self) -> &ModalPlateViscosity { &self.damping }
    /// Complete film memory on this equilibrium basis, when explicitly selected.
    /// Its film loss is excluded from `damping()` and must be attached exactly once.
    pub fn coating_relaxation(&self) -> Option<&super::coating::relaxation::ShellCoatingRelaxation> {
        self.coating_relaxation.as_ref()
    }
    /// Full causal substrate storage and loss in this same equilibrium basis.
    /// Its loss is excluded from `damping()`.
    pub fn substrate_relaxation(&self) -> Option<&super::constant_loss::relaxation::ShellSubstrateRelaxation> {
        self.substrate_relaxation.as_ref()
    }
    /// Work-conjugate interface rows in original source order and this mode basis.
    pub fn interface_columns(&self) -> &[Vec<f64>] { &self.interface_columns }
    /// Original interface viscosity as separate exact passive channels. Elastic
    /// stiffness is already in the monolithic pencil and is never applied twice.
    pub(crate) fn interface_damping_connections(&self, component: usize)
        -> Vec<crate::render::schedule::force::coupled::ModalConnection>
    {
        use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection};
        self.source.interfaces.iter().zip(&self.interface_columns)
            .filter(|(source, _)| source.law().damping_n_s_m != 0.0)
            .map(|(source, column)| ModalConnection {
                left:ModalAttachment { component, shapes:column.clone() },
                right:ModalAttachment { component, shapes:vec![0.0; column.len()] },
                stiffness_n_m:0.0, damping_n_s_m:source.law().damping_n_s_m,
                rest_extension_m:0.0,
            }).collect()
    }
    /// Compact monopole images; these have zero scalar loss, use C once.
    pub(crate) fn radiators(&self) -> &[CompactBody] { &self.radiators }
}
impl ResolvedShearShell {
    /// Compile one physical basis for arbitrary global force/moment footprints,
    /// signed volume change and every material-loss cross term. A force port's
    /// first three nodal values are force fractions; last three are moment/N
    /// lever arms [m], so the port remains displacement/force conjugate.
    ///
    /// Physical mode vectors are expanded by the ORIGINAL support-frame owner.
    /// Phi^T C Phi uses the elastic bubble transformation, not a second damped
    /// eigensolve. Retained modes and source bands are checked individually.
    /// No prestress or fluid loading is inferred. fs-modal's eigen call is
    /// synchronous; cancellation is checked immediately before and after it.
    pub fn reduce(&self, directors: &NodalDirectors, options: &ShellReductionOptions,
        force_ports: &[Vec<f64>], gate: &CancelGate) -> Result<MaterialShellImage, AcousticRealizeError>
    {
        poll(gate)?;
        let n=options.n_modes;
        let budget=options.projection_budget;
        let nodes=self.mesh.nodes.len(); let elements=self.mesh.tris.len();
        if !shell_mode_budget_admits(n,budget.max_modes)
            || !budget.reference_resistance_n_s_m.is_finite() || budget.reference_resistance_n_s_m<=0.0
            || !options.eigenvalue_window.0.is_finite() || !options.eigenvalue_window.1.is_finite()
            || options.eigenvalue_window.0<0.0 || options.eigenvalue_window.1<=options.eigenvalue_window.0
            || force_ports.is_empty() || force_ports.len()>4096
            || force_ports.iter().any(|p| p.len()!=6*nodes || p.iter().any(|v| !v.is_finite())
                || !p.iter().any(|v| *v!=0.0)) {
            return Err(invalid("shell reduction requires bounded modes, physical ports, window and resistance"));
        }
        let (exterior_facets,cavity_facets)=self.pressure_surfaces().map_or((elements,elements),
            |surfaces|(surfaces.exterior.len(),surfaces.cavity.len()));
        let work=self.loss_projection_terms(n)
            .and_then(|w|6usize.checked_mul(nodes)?.checked_mul(n)?.checked_mul(force_ports.len())?.checked_add(w))
            .and_then(|w|w.checked_add(4*n*n*n))
            .and_then(|w|w.checked_add(fs_plate::shell::pressure::projection_terms(nodes,exterior_facets,n)?))
            .and_then(|w|w.checked_add(fs_plate::shell::pressure::projection_terms(nodes,cavity_facets,n)?));
        if work.is_none_or(|w|w>budget.max_projection_terms) {
            return Err(invalid("complete shell condensation and physical projections exceed their declared allowance"));
        }
        let pressure=self.pressure_ports(gate)?;
        let model=assemble_mitc3p(&self.mesh,directors,&self.sections,&options.constraints,
            options.mass,options.assembly_budget,options.drilling_alpha,||gate.is_requested()).map_err(owner)?;
        let interface_laws: Vec<_> = self.interfaces.iter().map(|i| i.law().clone()).collect();
        let model = model.with_connections(&interface_laws, self.interface_budget,
            || gate.is_requested()).map_err(owner)?;
        // This is the shared spectral owner, not a private shell eigenalgorithm.
        let spectral=fs_modal::SliceOptions { mass_diagonal_equilibration:true, ..Default::default() };
        poll(gate)?;
        let report=modes_shell(model.pencil(),options.eigenvalue_window,&spectral).map_err(owner)?;
        poll(gate)?;
        if report.modes.len()<n {
            return Err(invalid("shell eigenvalue window contains fewer modes than requested"));
        }
        let modes=&report.modes[..n];
        let mut physical_modes=Vec::with_capacity(n);
        let mut intervals=Vec::with_capacity(n);
        let mut omega=Vec::with_capacity(n);
        for mode in modes {
            poll(gate)?;
            let w=fs_math::det::sqrt(mode.lambda);
            if !w.is_finite() || w<=0.0
                || self.regions.iter().any(|r|w<r.omega_band.0 || w>r.omega_band.1)
                || self.interfaces.iter().any(|r| {
                    let band = r.omega_band_rad_s(); w < band.0 || w > band.1
                }) {
                return Err(invalid("actual retained shell mode lies outside original material support"));
            }
            physical_modes.push(model.expand(&mode.phi).map_err(owner)?);
            omega.push(w); intervals.push(mode.interval);
        }
        let interface_columns = model.project_connections(&interface_laws, &physical_modes,
            self.interface_budget, || gate.is_requested()).map_err(owner)?;
        let mut projected=vec![vec![0.0;n];force_ports.len()];
        for (row,input) in projected.iter_mut().zip(force_ports) {
            poll(gate)?;
            // Independent work check in reduced support coordinates catches
            // accidental loss or double rotation of oblique constrained axes.
            let reduced=model.project_load(input).map_err(owner)?;
            for (j,physical) in physical_modes.iter().enumerate() {
                let value:f64=input.iter().zip(physical).map(|(a,b)|a*b).sum();
                let expected:f64=reduced.iter().zip(&modes[j].phi).map(|(a,b)|a*b).sum();
                let scale:f64=input.iter().zip(physical).map(|(a,b)|(a*b).abs()).sum();
                if !value.is_finite() || !scale.is_finite()
                    || (value-expected).abs()>1e-11*scale.max(1.0) {
                    return Err(invalid("shell physical force projection lost work conjugacy"));
                }
                row[j]=value;
            }
            if row.iter().all(|v|*v==0.0) {
                return Err(invalid("shell force footprint is detached from every retained mode"));
            }
        }
        let volume_shapes=pressure.exterior.project_modes(&physical_modes,
            budget.max_projection_terms,||gate.is_requested()).map_err(owner)?;
        let cavity_volume_shapes=pressure.cavity.project_modes(&physical_modes,
            budget.max_projection_terms,||gate.is_requested()).map_err(owner)?;
        let matrix=self.project_material_viscosity(directors,&physical_modes,&omega,
            options.assembly_budget,budget.max_projection_terms,gate)?;
        let damping=ModalPlateViscosity::factor(matrix,n,budget.reference_resistance_n_s_m,gate)?;
        let coating_relaxation=self.project_coating_relaxation(directors,&physical_modes,&omega,
            options.assembly_budget,budget.max_projection_terms,gate)?;
        let substrate_relaxation=self.project_substrate_relaxation(directors,&physical_modes,&omega,
            options.assembly_budget,budget.max_projection_terms,gate)?;
        let radiators=omega.iter().enumerate().map(|(i,&w)|
            CompactBody::from_mass_normalized_projection(w,volume_shapes[i],projected[0][i],pressure.exterior.area_m2()))
            .collect::<Result<Vec<_>,_>>()?;
        poll(gate)?;
        Ok(MaterialShellImage { source:self.clone(),directors:directors.clone(),physical_modes,force_ports:projected,
            volume_shapes,cavity_volume_shapes,intervals,radiators,damping,interface_columns,
            coating_relaxation,substrate_relaxation })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mode_budget_admits_exactly_up_to_the_ceiling() {
        // The former 64-mode limit still passes; the new ceiling passes; one more does not.
        assert!(shell_mode_budget_admits(1, 1));
        assert!(shell_mode_budget_admits(64, 64));
        assert!(shell_mode_budget_admits(65, 65));
        assert!(shell_mode_budget_admits(MAX_SHELL_MODES, MAX_SHELL_MODES));
        assert!(!shell_mode_budget_admits(MAX_SHELL_MODES + 1, MAX_SHELL_MODES + 1));
    }

    #[test]
    fn mode_budget_refuses_an_empty_request_or_one_above_its_own_declaration() {
        assert!(!shell_mode_budget_admits(0, MAX_SHELL_MODES));
        assert!(!shell_mode_budget_admits(3, 2));
        assert!(!shell_mode_budget_admits(1, 0));
        // A declared budget above the ceiling is refused even for a small request.
        assert!(!shell_mode_budget_admits(1, MAX_SHELL_MODES + 1));
    }
}
