//! Work-conjugate physical forcing and observations of the original wet frame.
//!
//! F_nodal=B*u, F_modal=Phi^T*B*u and y=B^T*Phi*v use the SAME real columns.
//! Ports may represent point forces, moments, or explicit distributed loads;
//! the caller owns units for general columns. Transverse member ports derive
//! their SI force/displacement map directly from the original structural trace.
//! No normalization, absolute-value replacement or fitted pressure is introduced.

use super::MaterialFluidFrame;
use super::harmonic::{FrameHarmonicError, FrameHarmonicResponse, PreparedFrameHarmonic};
use fs_exec::CancelGate;
use fs_math::c64::C64;

#[derive(Clone, Copy, Debug)]
pub struct FrameHarmonicPortBudget {
    pub max_ports: usize,
    /// ports * structural_modes * original_nodal_dofs, before projection.
    pub max_projection_terms: usize,
    /// Entries in the complete ports-by-ports mobility matrix.
    pub max_mobility_entries: usize,
    pub power_absolute_tolerance_w: f64,
    pub power_relative_tolerance: f64,
}
#[derive(Debug)]
pub enum FrameHarmonicPortError {
    Invalid(&'static str),
    Harmonic(FrameHarmonicError),
    Structural(fs_solid::SolidError),
    PowerMapping { residual_w: f64, tolerance_w: f64 },
    Cancelled,
}
impl core::fmt::Display for FrameHarmonicPortError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f,"frame harmonic ports: {self:?}")
    }
}
impl std::error::Error for FrameHarmonicPortError {}
impl From<FrameHarmonicError> for FrameHarmonicPortError {
    fn from(error: FrameHarmonicError) -> Self {
        match error { FrameHarmonicError::Cancelled => Self::Cancelled, other => Self::Harmonic(other) }
    }
}
type Result<T> = core::result::Result<T,FrameHarmonicPortError>;
fn invalid(message: &'static str) -> FrameHarmonicPortError { FrameHarmonicPortError::Invalid(message) }
fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() { Err(FrameHarmonicPortError::Cancelled) } else { Ok(()) }
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() { Ok(x) } else { Err(invalid("unrepresentable physical port arithmetic")) }
}
fn complex(x: C64) -> Result<C64> { finite(x.re)?; finite(x.im)?; Ok(x) }
fn projection(column: &[f64], values: &[C64]) -> Result<C64> {
    column.iter().zip(values).try_fold(C64::ZERO,|s,(&b,&x)|complex(s+x.scale(b)))
}

/// Source-borrowing physical columns, projected once for an entire sweep.
/// A factor prepared from any OTHER model refuses, even if dimensions match.
/// Repeated, signed, dependent and zero columns remain meaningful; they are
/// never normalized or made artificially independent.
pub struct FrameHarmonicPorts<'m,'f,F> {
    model: &'m MaterialFluidFrame<'f,F>,
    nodal_columns: Vec<Vec<f64>>,
    modal_columns: Vec<Vec<f64>>,
    budget: FrameHarmonicPortBudget,
}

#[derive(Clone, Debug)]
pub struct FrameHarmonicPortResponse {
    pub mechanical: FrameHarmonicResponse,
    /// Peak phasors in original node-major (x,y,rotation) order. Displacement is
    /// the harmonic increment, NOT the frozen equilibrium/reference offset.
    pub nodal_displacement: Vec<C64>,
    pub nodal_velocity: Vec<C64>,
    pub nodal_force: Vec<C64>,
    pub port_displacement: Vec<C64>,
    pub port_velocity: Vec<C64>,
    /// Signed mean work per drive. A drive may ABSORB work from other drives.
    pub port_power_w: Vec<f64>,
    pub nodal_power_w: f64,
    pub mapping_power_residual_w: f64,
}

/// Row-major velocity-per-effort matrix with peak-phasor exp(+i omega t).
/// No symmetry averaging, gain fitting, modal normalization or passivity clamp.
/// This steady-response matrix is not itself a causal time-domain realization.
#[derive(Clone, Debug)]
pub struct FrameMobilityMatrix {
    pub omega_rad_s: f64,
    pub port_count: usize,
    pub entries: Vec<C64>,
    pub maximum_equation_relative_residual: f64,
    pub maximum_power_residual_w: f64,
}

fn admit(modes: usize, dofs: usize, ports: usize, budget: FrameHarmonicPortBudget) -> Result<()> {
    if ports == 0 || ports > budget.max_ports
        || ports.checked_mul(modes).and_then(|n|n.checked_mul(dofs))
            .is_none_or(|n| n>budget.max_projection_terms)
        || ports.checked_mul(ports).is_none_or(|n| n>budget.max_mobility_entries)
        || !budget.power_absolute_tolerance_w.is_finite() || budget.power_absolute_tolerance_w<=0.0
        || !budget.power_relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&budget.power_relative_tolerance) {
        return Err(invalid("complete physical port map and finite work budget required"));
    }
    Ok(())
}

impl<'f,F> MaterialFluidFrame<'f,F> {
    /// Bind original node-major force/moment columns to their conjugate motion.
    /// Each column has exactly one coefficient per original nodal DOF. Units
    /// must make column*effort a force [N] or moment [N m] at the relevant row.
    /// Constraint reactions do no work because their retained displacement rows
    /// are zero. No force or observation is applied to private fluid memories.
    pub fn harmonic_ports(&self, columns: &[Vec<f64>], budget: FrameHarmonicPortBudget,
        gate: &CancelGate) -> Result<FrameHarmonicPorts<'_,'f,F>> {
        poll(gate)?;
        let m=self.basis.modes().len();
        let dofs=self.basis.reference_displacement().len();
        admit(m,dofs,columns.len(),budget)?;
        for column in columns {
            poll(gate)?;
            if column.len()!=dofs || column.iter().any(|x|!x.is_finite()) {
                return Err(invalid("physical port column needs every finite nodal force/moment coefficient"));
            }
        }
        let mut projected=Vec::with_capacity(columns.len());
        for column in columns {
            let mut modal=Vec::with_capacity(m);
            for mode in self.basis.modes() {
                poll(gate)?;
                let mut value=0.0;
                for (&b,&shape) in column.iter().zip(&mode.displacement) { value=finite(value+b*shape)?; }
                modal.push(value);
            }
            projected.push(modal);
        }
        poll(gate)?;
        Ok(FrameHarmonicPorts { model:self,nodal_columns:columns.to_vec(),modal_columns:projected,budget })
    }

    /// Point-force ports at (original member index, reference xi). The original
    /// cubic transverse trace supplies BOTH nodal forces/moments and observed
    /// displacement; no separate point geometry or endpoint-only map is invented.
    /// Positive effort is along each member's original reference normal.
    pub fn harmonic_transverse_ports(&self, points: &[(usize,f64)],
        budget: FrameHarmonicPortBudget, gate: &CancelGate) -> Result<FrameHarmonicPorts<'_,'f,F>> {
        poll(gate)?;
        let dofs=self.basis.reference_displacement().len();
        admit(self.basis.modes().len(),dofs,points.len(),budget)?;
        let mut columns=Vec::with_capacity(points.len());
        for &(element,xi) in points {
            poll(gate)?;
            let (indices,weights,_)=self.frame.reference_transverse_trace(element,xi)
                .map_err(FrameHarmonicPortError::Structural)?;
            let mut column=vec![0.0;dofs];
            for (&index,&weight) in indices.iter().zip(&weights) { column[index]=finite(column[index]+weight)?; }
            columns.push(column);
        }
        self.harmonic_ports(&columns,budget,gate)
    }
}

impl<'m,'f,F> FrameHarmonicPorts<'m,'f,F> {
    pub const fn model(&self) -> &'m MaterialFluidFrame<'f,F> { self.model }
    pub fn nodal_columns(&self) -> &[Vec<f64>] { &self.nodal_columns }
    pub fn modal_columns(&self) -> &[Vec<f64>] { &self.modal_columns }
    pub fn port_count(&self) -> usize { self.nodal_columns.len() }
    fn admit_factor(&self, prepared: &PreparedFrameHarmonic<'_,'_,F>) -> Result<()> {
        if !core::ptr::eq(self.model,prepared.model()) {
            return Err(invalid("physical port map and frequency factor require the SAME original model"));
        }
        Ok(())
    }

    /// Combine complex drive amplitudes BEFORE the physical solve and loss
    /// calculation. Summing separate unit-drive heats would omit interference.
    /// Reports and native triangular solves allocate; source preparation and
    /// complex factorization are reused. No accepted history is advanced.
    pub fn solve(&self, prepared: &PreparedFrameHarmonic<'_,'_,F>, efforts: &[C64],
        gate: &CancelGate) -> Result<FrameHarmonicPortResponse> {
        poll(gate)?;
        self.admit_factor(prepared)?;
        if efforts.len()!=self.port_count() { return Err(invalid("one complex amplitude per physical drive required")); }
        let m=self.model.basis.modes().len();
        let dofs=self.model.basis.reference_displacement().len();
        let mut force=vec![C64::ZERO;m];
        let mut nodal_force=vec![C64::ZERO;dofs];
        for (port,&effort) in efforts.iter().enumerate() {
            poll(gate)?; complex(effort)?;
            for (value,&b) in force.iter_mut().zip(&self.modal_columns[port]) { *value=complex(*value+effort.scale(b))?; }
            for (value,&b) in nodal_force.iter_mut().zip(&self.nodal_columns[port]) { *value=complex(*value+effort.scale(b))?; }
        }
        let mechanical=prepared.solve_modal(&force,gate)?;
        let mut displacement=vec![C64::ZERO;dofs];
        let mut velocity=vec![C64::ZERO;dofs];
        for (mode,source) in self.model.basis.modes().iter().enumerate() {
            poll(gate)?;
            for dof in 0..dofs {
                displacement[dof]=complex(displacement[dof]+mechanical.modal_displacement[mode].scale(source.displacement[dof]))?;
                velocity[dof]=complex(velocity[dof]+mechanical.modal_velocity[mode].scale(source.displacement[dof]))?;
            }
        }
        let mut port_displacement=Vec::with_capacity(efforts.len());
        let mut port_velocity=Vec::with_capacity(efforts.len());
        let mut port_power=Vec::with_capacity(efforts.len());
        for (column,&effort) in self.nodal_columns.iter().zip(efforts) {
            poll(gate)?;
            port_displacement.push(projection(column,&displacement)?);
            let speed=projection(column,&velocity)?;
            port_velocity.push(speed);
            port_power.push(finite(0.5*complex(effort*speed.conj())?.re)?);
        }
        let nodal_power=nodal_force.iter().zip(&velocity).try_fold(0.0,|s,(&f,&v)|
            finite(s+0.5*complex(f*v.conj())?.re))?;
        let total_port=port_power.iter().try_fold(0.0,|s,&p|finite(s+p))?;
        let defects=[finite(nodal_power-mechanical.supplied_power_w)?,finite(total_port-nodal_power)?];
        let residual=if defects[0].abs()>=defects[1].abs() {defects[0]} else {defects[1]};
        let scale=nodal_power.abs().max(total_port.abs()).max(mechanical.supplied_power_w.abs());
        let tolerance=finite(self.budget.power_absolute_tolerance_w+self.budget.power_relative_tolerance*scale)?;
        if residual.abs()>tolerance { return Err(FrameHarmonicPortError::PowerMapping {residual_w:residual,tolerance_w:tolerance}); }
        poll(gate)?;
        Ok(FrameHarmonicPortResponse { mechanical,nodal_displacement:displacement,nodal_velocity:velocity,
            nodal_force,port_displacement,port_velocity,port_power_w:port_power,
            nodal_power_w:nodal_power,mapping_power_residual_w:residual })
    }

    /// Complete Y_ij = velocity at port i / unit effort at port j. The supplied
    /// factor is reused for EVERY column; no per-port eigensolve, quadrature or
    /// factorization occurs. All source losses and mapping gates run per column.
    /// Rank-deficient port layouts are valid and are not artificially repaired.
    /// For coherent drives use `solve`, not a sum of independent column heats.
    pub fn mobility(&self, prepared: &PreparedFrameHarmonic<'_,'_,F>, gate: &CancelGate)
        -> Result<FrameMobilityMatrix> {
        poll(gate)?; self.admit_factor(prepared)?;
        let count=self.port_count();
        let mut entries=vec![C64::ZERO;count*count];
        let mut efforts=vec![C64::ZERO;count];
        let mut maximum_equation_relative_residual=0.0_f64;
        let mut maximum_power_residual_w=0.0_f64;
        for column in 0..count {
            poll(gate)?;efforts[column]=C64::new(1.0,0.0);
            let response=self.solve(prepared,&efforts,gate)?;
            efforts[column]=C64::ZERO;
            for row in 0..count {entries[row*count+column]=response.port_velocity[row];}
            maximum_equation_relative_residual=maximum_equation_relative_residual.max(response.mechanical.equation_relative_residual);
            maximum_power_residual_w=maximum_power_residual_w.max(response.mechanical.power_residual_w.abs())
                .max(response.mapping_power_residual_w.abs());
        }
        poll(gate)?;
        Ok(FrameMobilityMatrix {omega_rad_s:prepared.omega_rad_s(),port_count:count,entries,
            maximum_equation_relative_residual,maximum_power_residual_w})
    }
}
