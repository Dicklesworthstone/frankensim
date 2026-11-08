//! Coherent finite-time heat and Lorentz forcing on the original 3D AC fields.
//!
//! Real fields are `Re(sum field_j exp(-i omega_j t))`. Products retain both
//! difference-frequency beats and sum-frequency carriers. Each line keeps its
//! resolved skin/proximity distribution; averaging each line independently
//! would lose their interference. Equal frequencies combine before products.
//! Exact tetrahedral P1 mass weights and the existing analytic temporal Gram
//! kernel perform the integration without temporal sampling or another solve.
//!
//! Every response must belong to the SAME unchanged source and prescribed
//! steady velocity. The source's scalar/oriented tensor conductivity therefore
//! remains the constitutive owner. An optional admitted static magnetic bias
//! contributes its carrier-frequency force on stationary AC conduction current.
//! DC current forces/heat, switch-on transients, nonlinear/hysteretic mixing,
//! moving geometry and mechanical feedback remain outside this owner. Heat must be
//! deposited once, without adding the original cycle means again. No new
//! continuum error or complete source/storage-balance certificate is supplied.

use super::super::mechanical::{
    SpatialLorentzLoadError, SpatialLorentzLoadLimits, SpatialLorentzLoads,
    load::{CellLoad, integrate},
};
use super::{
    SpatialHarmonicInductionResponse, SpatialInduction,
    loads::{admit_static_bias, admit_static_bias_cell},
};
use crate::spatial::vector_potential::SpatialVectorField;
use fs_fft::harmonic_window::{
    HarmonicProduct, HarmonicWindow, HarmonicWindowProfile, WindowError,
};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_math::{c64::C64, dd::Dd};

/// Hard ceiling; the caller's pair and workspace budgets also apply.
pub const MAX_EXPOSURE_LINES: usize = 256;

/// Additional work/storage admission, excluding the borrowed native fields.
#[derive(Clone, Copy, Debug)]
pub struct SpatialHarmonicInductionExposureLimits {
    pub maximum_lines: usize,
    /// Conservative scalar temporal contractions per observation: 39 for each
    /// original conductor and ordered pair of authored lines. Equal-frequency
    /// grouping may reduce the actual work, but never bypasses this admission.
    /// A DC-biased observation additionally needs 12 contractions for each
    /// original conductor and authored line, within this same bound.
    pub maximum_pair_products: usize,
    /// Prepared bank, velocity, grouping, cached pair kernels and cell output.
    /// Allocator metadata and constant-size stack scratch are excluded.
    /// A DC-biased observation also caches one moment kernel per authored line;
    /// the original static field is borrowed without a field-sized copy.
    /// Optional nodal loads have their own `SpatialLorentzLoadLimits`.
    pub maximum_workspace_bytes: usize,
}

/// No partial observation or change of accepted field state on refusal.
#[derive(Debug)]
pub enum SpatialHarmonicInductionExposureError {
    Invalid(&'static str),
    Cancelled,
    Moment(WindowError),
    Loads(SpatialLorentzLoadError),
}
impl core::fmt::Display for SpatialHarmonicInductionExposureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "3D harmonic induction exposure: {self:?}")
    }
}
impl std::error::Error for SpatialHarmonicInductionExposureError {}
impl From<WindowError> for SpatialHarmonicInductionExposureError {
    fn from(error: WindowError) -> Self {
        Self::Moment(error)
    }
}
impl From<SpatialLorentzLoadError> for SpatialHarmonicInductionExposureError {
    fn from(error: SpatialLorentzLoadError) -> Self {
        match error {
            SpatialLorentzLoadError::Cancelled => Self::Cancelled,
            error => Self::Loads(error),
        }
    }
}
type Result<T> = core::result::Result<T, SpatialHarmonicInductionExposureError>;
type Phasor = [Dd; 2];
type Vector = [Phasor; 3];

/// One original conductor's AC heat and finite-interval mean force density.
pub struct SpatialHarmonicInductionCellExposure<'s, C: SpatialConductorMaterial> {
    pub cell: usize,
    pub source: &'s C,
    pub volume_m3: f64,
    pub mean_ohmic_power_w: f64,
    pub ohmic_heat_j: f64,
    /// P1 values at the original four tetrahedron vertices, in field coordinates.
    pub mean_lorentz_force_density_n_m3: [[f64; 3]; 4],
    /// Signed field-to-prescribed-motion power, separate from thermal deposition.
    pub mean_mechanical_power_w: f64,
}

/// An unweighted physical interval on the common `exp(-i omega t)` source clock.
pub struct SpatialHarmonicInductionExposureResult<
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    source: &'s SpatialInduction<'m, C>,
    cells: Vec<SpatialHarmonicInductionCellExposure<'s, C>>,
    pub window: HarmonicWindow,
    pub mean_ohmic_power_w: f64,
    pub ohmic_heat_j: f64,
    pub mean_mechanical_power_w: f64,
    pub mechanical_work_j: f64,
}

struct LineCell {
    electric: [Vector; 4],
    current: [Vector; 4],
    magnetic: Vector,
    electric_sum: Vector,
    current_sum: Vector,
}
impl LineCell {
    fn zero() -> Self {
        Self {
            electric: [[[Dd::ZERO; 2]; 3]; 4],
            current: [[[Dd::ZERO; 2]; 3]; 4],
            magnetic: [[Dd::ZERO; 2]; 3],
            electric_sum: [[Dd::ZERO; 2]; 3],
            current_sum: [[Dd::ZERO; 2]; 3],
        }
    }
}
struct Line {
    omega_rad_s: f64,
    cells: Vec<LineCell>,
}

/// Reusable original fields, with O(F² * conducting_cells) observation work.
/// Pair trigonometry is evaluated once per window, never per spatial point.
pub struct SpatialHarmonicInductionExposure<
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    source: &'s SpatialInduction<'m, C>,
    velocity_m_s: Vec<[f64; 3]>,
    lines: Vec<Line>,
    authored_lines: usize,
    limits: SpatialHarmonicInductionExposureLimits,
    base_products: usize,
    base_workspace_bytes: usize,
}

fn invalid(message: &'static str) -> SpatialHarmonicInductionExposureError {
    SpatialHarmonicInductionExposureError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(SpatialHarmonicInductionExposureError::Cancelled)
    } else {
        Ok(())
    }
}
fn real(value: Dd) -> Result<f64> {
    let rounded = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && rounded.is_finite() {
        Ok(rounded)
    } else {
        Err(invalid("harmonic exposure arithmetic is not representable"))
    }
}
fn nonnegative(value: Dd) -> Result<f64> {
    let value = real(value)?;
    if value < 0.0 {
        Err(invalid("negative coherent conductor heat"))
    } else {
        Ok(value)
    }
}
fn add(sum: &mut Phasor, value: C64) {
    sum[0] = sum[0] + Dd::from_f64(value.re);
    sum[1] = sum[1] + Dd::from_f64(value.im);
}

impl<'m, C: SpatialConductorMaterial> SpatialInduction<'m, C> {
    /// Prepare finite-window observations of already accepted native AC solves.
    ///
    /// All lines share this exact source, material state, geometry and prescribed
    /// velocity. Frequencies and peak phases stay those of the original solves;
    /// a common source clock is the caller's physical declaration. Exact equal
    /// frequencies sum in DD before squaring, including their low components.
    /// Distinct frequencies are never merged by a tolerance. Preparation and
    /// observation do not mutate histories, advance clocks or solve new fields.
    pub fn prepare_harmonic_exposure<'s>(
        &'s self,
        responses: &[&SpatialHarmonicInductionResponse<'_, 'm, C>],
        limits: SpatialHarmonicInductionExposureLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionExposure<'s, 'm, C>> {
        poll(cancelled)?;
        let count = responses.len();
        if count == 0 || count > MAX_EXPOSURE_LINES || count > limits.maximum_lines {
            return Err(invalid(
                "nonempty bounded original AC response list required",
            ));
        }
        let cells = self.prepared.cells.len();
        let vertices = self.problem().mesh().vertex_count;
        let pairs = count.checked_mul(count);
        let products = pairs
            .and_then(|n| n.checked_mul(cells))
            .and_then(|n| n.checked_mul(39))
            .ok_or_else(|| invalid("3D exposure work size overflow"))?;
        let bytes = cells
            .checked_mul(core::mem::size_of::<LineCell>())
            .and_then(|n| n.checked_add(core::mem::size_of::<Line>()))
            .and_then(|n| n.checked_mul(count))
            .and_then(|n| {
                pairs
                    .and_then(|v| v.checked_mul(core::mem::size_of::<HarmonicProduct>()))
                    .and_then(|v| n.checked_add(v))
            })
            .and_then(|n| {
                vertices
                    .checked_mul(core::mem::size_of::<[f64; 3]>())
                    .and_then(|v| n.checked_add(v))
            })
            .and_then(|n| {
                cells
                    .checked_mul(core::mem::size_of::<
                        SpatialHarmonicInductionCellExposure<'_, C>,
                    >())
                    .and_then(|v| n.checked_add(v))
            })
            .and_then(|n| {
                count
                    .checked_mul(core::mem::size_of::<usize>())
                    .and_then(|v| n.checked_add(v))
            })
            .and_then(|n| {
                n.checked_add(core::mem::size_of::<
                    SpatialHarmonicInductionExposure<'_, '_, C>,
                >())
            })
            .ok_or_else(|| invalid("3D exposure workspace size overflow"))?;
        if products > limits.maximum_pair_products || bytes > limits.maximum_workspace_bytes {
            return Err(invalid("3D exposure exceeds pair or workspace budget"));
        }
        let first = responses[0];
        for response in responses {
            poll(cancelled)?;
            if !core::ptr::eq(response.source, self) || response.velocity_m_s != first.velocity_m_s
            {
                return Err(invalid(
                    "3D exposure requires the same original source and steady velocity",
                ));
            }
        }
        let mut order: Vec<_> = (0..count).collect();
        order.sort_unstable_by(|&a, &b| {
            responses[a]
                .omega_rad_s
                .total_cmp(&responses[b].omega_rad_s)
                .then_with(|| a.cmp(&b))
        });
        let mut lines = Vec::with_capacity(count);
        let mut start = 0;
        while start < count {
            poll(cancelled)?;
            let omega_rad_s = responses[order[start]].omega_rad_s;
            let mut end = start + 1;
            while end < count && responses[order[end]].omega_rad_s == omega_rad_s {
                end += 1;
            }
            let mut fields = Vec::with_capacity(cells);
            for index in 0..cells {
                let mut field = LineCell::zero();
                for &response_index in &order[start..end] {
                    poll(cancelled)?;
                    let response = responses[response_index];
                    let record = &response.cells[index];
                    for axis in 0..3 {
                        add(
                            &mut field.magnetic[axis],
                            response.flux_density_t[record.cell][axis],
                        );
                        for vertex in 0..4 {
                            add(
                                &mut field.electric[vertex][axis],
                                record.conductor_electric_field_v_m[vertex][axis],
                            );
                            add(
                                &mut field.current[vertex][axis],
                                record.current_density_a_m2[vertex][axis],
                            );
                        }
                    }
                }
                for axis in 0..3 {
                    for part in 0..2 {
                        real(field.magnetic[axis][part])?;
                        for vertex in 0..4 {
                            real(field.electric[vertex][axis][part])?;
                            real(field.current[vertex][axis][part])?;
                            field.electric_sum[axis][part] =
                                field.electric_sum[axis][part] + field.electric[vertex][axis][part];
                            field.current_sum[axis][part] =
                                field.current_sum[axis][part] + field.current[vertex][axis][part];
                        }
                        real(field.electric_sum[axis][part])?;
                        real(field.current_sum[axis][part])?;
                    }
                }
                fields.push(field);
            }
            lines.push(Line {
                omega_rad_s,
                cells: fields,
            });
            start = end;
        }
        poll(cancelled)?;
        Ok(SpatialHarmonicInductionExposure {
            source: self,
            velocity_m_s: first.velocity_m_s.clone(),
            lines,
            authored_lines: count,
            limits,
            base_products: products,
            base_workspace_bytes: bytes,
        })
    }
}

impl<'s, 'm, C: SpatialConductorMaterial> SpatialHarmonicInductionExposure<'s, 'm, C> {
    pub fn source(&self) -> &'s SpatialInduction<'m, C> {
        self.source
    }
    pub fn distinct_frequencies(&self) -> usize {
        self.lines.len()
    }
    /// Integrate real J·E_conductor and J×B over one physical time interval.
    ///
    /// Only uniform windows give physical deposited heat and mechanical work;
    /// shaped detector windows are refused. All original spatial cross terms
    /// and both temporal sidebands survive. Small negative heat caused by
    /// unresolvable cancellation is refused, not clipped into a physical result.
    /// Analytic time integration removes sampling error but does not remove the
    /// original field-solve, mesh, material or floating-point errors.
    pub fn observe(
        &self,
        window: HarmonicWindow,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionExposureResult<'s, 'm, C>> {
        self.observe_inner(window, None, cancelled)
    }

    /// Integrate the AC fields with the additional force `J_ac(t) cross B_dc`.
    ///
    /// The admitted native static field must share the original magnetic
    /// problem and prescribed-edge support. All AC responses must have zero
    /// prescribed velocity: the static magnetic owner supplies no DC motional
    /// current. All original conducting cells must also have zero impressed
    /// DC current, preserving the native static-bias load admission. Static
    /// sources in insulating cells may supply the field. No transient
    /// induction state is relabeled as a steady bias.
    ///
    /// Each line's mean current uses the shared analytic zero-carrier moment,
    /// retaining phase and finite-window carrier impulse without sampling.
    /// Grouped DD currents combine before multiplication by the static field.
    /// The original preparation limits must also admit 12*C*F additional
    /// temporal contractions and F additional cached kernels, using AUTHORED
    /// line count F even when equal-frequency lines coalesced. Refusal precedes
    /// the extra workspace allocation; the native DC field is only borrowed.
    ///
    /// AC Joule heat is unchanged. This result concerns the original resolved
    /// AC conduction current only. Impressed DC source forces/heat and their
    /// AC interaction are outside this observation, and DC-driven conducting
    /// cells refuse rather than omitting their additional force. DC terminal
    /// currents, source/storage closure and mechanical feedback also require
    /// their own coupled owners.
    pub fn observe_with_static_bias(
        &self,
        bias: &SpatialVectorField<'_>,
        window: HarmonicWindow,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionExposureResult<'s, 'm, C>> {
        poll(cancelled)?;
        let products = self
            .authored_lines
            .checked_mul(self.source.prepared.cells.len())
            .and_then(|n| n.checked_mul(12))
            .and_then(|n| n.checked_add(self.base_products));
        let bytes = self
            .authored_lines
            .checked_mul(core::mem::size_of::<HarmonicProduct>())
            .and_then(|n| n.checked_add(self.base_workspace_bytes));
        if products.is_none_or(|n| n > self.limits.maximum_pair_products)
            || bytes.is_none_or(|n| n > self.limits.maximum_workspace_bytes)
        {
            return Err(invalid(
                "DC-biased exposure exceeds its original work or workspace budget",
            ));
        }
        admit_static_bias(self.source, &self.velocity_m_s, bias, cancelled)?;
        for cell in &self.source.prepared.cells {
            poll(cancelled)?;
            admit_static_bias_cell(bias, cell.cell)?;
        }
        self.observe_inner(window, Some(bias), cancelled)
    }

    fn observe_inner(
        &self,
        window: HarmonicWindow,
        bias: Option<&SpatialVectorField<'_>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicInductionExposureResult<'s, 'm, C>> {
        poll(cancelled)?;
        if window.profile() != HarmonicWindowProfile::Uniform {
            return Err(invalid(
                "physical heat and work require a uniform time interval",
            ));
        }
        let count = self.lines.len();
        let mut kernels = Vec::with_capacity(count * count);
        for a in &self.lines {
            for b in &self.lines {
                poll(cancelled)?;
                kernels.push(window.product_kernel(a.omega_rad_s, b.omega_rad_s)?);
            }
        }
        let mut bias_kernels = Vec::with_capacity(if bias.is_some() { count } else { 0 });
        if bias.is_some() {
            for line in &self.lines {
                poll(cancelled)?;
                bias_kernels.push(window.product_kernel(line.omega_rad_s, 0.0)?);
            }
        }
        let source: &'s SpatialInduction<'m, C> = self.source;
        let mut cells = Vec::with_capacity(source.prepared.cells.len());
        let mut heat = Dd::ZERO;
        let mut mechanical = Dd::ZERO;
        let duration = Dd::from_f64(window.duration_s());
        for (index, original) in source.prepared.cells.iter().enumerate() {
            poll(cancelled)?;
            let mut product = Dd::ZERO;
            let mut force = [[Dd::ZERO; 3]; 4];
            for (i, a) in self.lines.iter().enumerate() {
                for (j, b) in self.lines.iter().enumerate() {
                    poll(cancelled)?;
                    let (a, b, kernel) = (&a.cells[index], &b.cells[index], kernels[i * count + j]);
                    for axis in 0..3 {
                        product = product
                            + kernel.mean_real_product_parts_dd(
                                a.current_sum[axis],
                                b.electric_sum[axis],
                            )?;
                        let (next, last) = ((axis + 1) % 3, (axis + 2) % 3);
                        for vertex in 0..4 {
                            product = product
                                + kernel.mean_real_product_parts_dd(
                                    a.current[vertex][axis],
                                    b.electric[vertex][axis],
                                )?;
                            force[vertex][axis] = force[vertex][axis]
                                + kernel.mean_real_product_parts_dd(
                                    a.current[vertex][next],
                                    b.magnetic[last],
                                )?
                                - kernel.mean_real_product_parts_dd(
                                    a.current[vertex][last],
                                    b.magnetic[next],
                                )?;
                        }
                    }
                }
            }
            if let Some(bias) = bias {
                let magnetic = bias.flux_density_t()[original.cell].map(Dd::from_f64);
                let mut mean_current = [[Dd::ZERO; 3]; 4];
                for (line, kernel) in self.lines.iter().zip(&bias_kernels) {
                    poll(cancelled)?;
                    for vertex in 0..4 {
                        for axis in 0..3 {
                            mean_current[vertex][axis] = mean_current[vertex][axis]
                                + kernel.mean_real_product_parts_dd(
                                    line.cells[index].current[vertex][axis],
                                    [Dd::from_f64(1.0), Dd::ZERO],
                                )?;
                        }
                    }
                }
                for vertex in 0..4 {
                    for axis in 0..3 {
                        let (next, last) = ((axis + 1) % 3, (axis + 2) % 3);
                        force[vertex][axis] = force[vertex][axis]
                            + mean_current[vertex][next] * magnetic[last]
                            - mean_current[vertex][last] * magnetic[next];
                    }
                }
            }
            // Integral N_a N_b = V*(1+delta_ab)/20, without spatial lumping.
            let mass = Dd::from_f64(original.volume) / Dd::from_f64(20.0);
            let watts = mass * product;
            let mean_ohmic_power_w = nonnegative(watts)?;
            let tet = source.problem().mesh().tets[original.cell];
            let mut work = Dd::ZERO;
            let mut force_density = [[0.0; 3]; 4];
            for axis in 0..3 {
                let mut sum_force = Dd::ZERO;
                let mut sum_velocity = Dd::ZERO;
                for vertex in 0..4 {
                    let velocity = Dd::from_f64(self.velocity_m_s[tet[vertex] as usize][axis]);
                    force_density[vertex][axis] = real(force[vertex][axis])?;
                    sum_force = sum_force + force[vertex][axis];
                    sum_velocity = sum_velocity + velocity;
                    work = work + force[vertex][axis] * velocity;
                }
                work = work + sum_force * sum_velocity;
            }
            work = mass * work;
            heat = heat + watts;
            mechanical = mechanical + work;
            cells.push(SpatialHarmonicInductionCellExposure {
                cell: original.cell,
                source: source.conductors()[original.cell]
                    .as_ref()
                    .ok_or_else(|| invalid("missing original conductor source"))?,
                volume_m3: original.volume,
                mean_ohmic_power_w,
                ohmic_heat_j: nonnegative(duration * watts)?,
                mean_lorentz_force_density_n_m3: force_density,
                mean_mechanical_power_w: real(work)?,
            });
        }
        poll(cancelled)?;
        Ok(SpatialHarmonicInductionExposureResult {
            source,
            cells,
            window,
            mean_ohmic_power_w: nonnegative(heat)?,
            ohmic_heat_j: nonnegative(duration * heat)?,
            mean_mechanical_power_w: real(mechanical)?,
            mechanical_work_j: real(duration * mechanical)?,
        })
    }
}

impl<'s, 'm, C: SpatialConductorMaterial> SpatialHarmonicInductionExposureResult<'s, 'm, C> {
    pub fn source(&self) -> &'s SpatialInduction<'m, C> {
        self.source
    }
    /// Immutable original-order cells, with their exact source state retained.
    pub fn cells(&self) -> &[SpatialHarmonicInductionCellExposure<'s, C>] {
        &self.cells
    }
    /// Consistent nodal mean loads for distinct selected original conductors.
    /// The returned native wrench and P1 virtual-power APIs use the original
    /// coordinates. Multiplying a mean force/wrench by this window's duration
    /// gives its impulse. This transfer supplies one-way forcing, not a closed
    /// mechanical evolution, and excludes overlapping Maxwell surface loads.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialLorentzLoads> {
        poll(cancelled)?;
        if selected_cells.is_empty() || selected_cells.len() > limits.max_cells {
            return Err(invalid("nonempty bounded conductor selection required"));
        }
        let mut selected = Vec::with_capacity(selected_cells.len());
        for &cell in selected_cells {
            poll(cancelled)?;
            let index = self
                .cells
                .binary_search_by_key(&cell, |record| record.cell)
                .map_err(|_| invalid("exposure load selection requires original conductors"))?;
            selected.push(&self.cells[index]);
        }
        selected.sort_unstable_by_key(|record| record.cell);
        if selected.windows(2).any(|pair| pair[0].cell == pair[1].cell) {
            return Err(invalid("repeated exposure load cell"));
        }
        let problem = self.source.problem();
        Ok(integrate(
            problem.positions_m(),
            selected.iter().map(|record| CellLoad {
                vertices: problem.mesh().tets[record.cell],
                volume_m3: record.volume_m3,
                force_density_n_m3: &record.mean_lorentz_force_density_n_m3,
            }),
            limits,
            cancelled,
        )?)
    }
}
