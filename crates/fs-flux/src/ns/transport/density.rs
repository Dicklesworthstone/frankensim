//! Conservative P0 density and specific-inventory transport on native BDM faces.
//!
//! A prescribed divergence-free volume velocity transports mass, not a fixed
//! cell inventory. Each donor's density multiplies its OWN directional volume
//! rate, including simultaneous opposing streams. The resulting face mass
//! transfers also carry specific quantities, e.g. [u_x, u_y, h]. This advances
//! cell momentum and enthalpy inventories with exactly the same mass transfer;
//! it does not solve the BDM momentum equations, pressure, an EOS, or heat flow.
//! In particular it is not a changing-density Navier–Stokes time integrator.
//!
//! Synchronous donor-cell Euler is first order and numerically mixes material.
//! Positive remaining donor mass requires an outgoing-volume Courant number
//! strictly below one. No density clipping, mass repair, or interface sharpening
//! is performed. Fixed volumes and divergence-free prescribed flow exclude
//! compressibility, phase-volume sources and moving-mesh geometric conservation.

use super::{EdgeVolumeTransport, FluxTransportError, dimensional_rate, edge_volume_rate};
use crate::{FluxSystem, TriMesh};

/// Cell-average density and N quantities per unit mass. Components may have
/// different units: velocity gives momentum, specific enthalpy gives enthalpy,
/// and mass fractions give component masses. Their reference convention and
/// constitutive admissibility remain with the caller. N=0 transports only mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityTransportCell<const N: usize> {
    /// Positive mass density [kg/m³] in the original fixed cell volume.
    pub density_kg_m3: f64,
    /// Finite specific values, each with its caller-declared per-kg convention.
    pub specific: [f64; N],
}

impl<const N: usize> Default for DensityTransportCell<N> {
    fn default() -> Self {
        Self {
            density_kg_m3: 0.0,
            specific: [0.0; N],
        }
    }
}

/// Time-integrated transfer through one original edge. The global normal is
/// outward from `edge.tris.0`; on a boundary the reverse stream is inflow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EdgeDensityTransfer<const N: usize> {
    /// Nonnegative directional mass amounts [kg], never collapsed to net flow.
    pub mass_kg: [f64; 2],
    /// Extensive carried amounts in the same [along, against] direction order.
    /// A component can be signed, e.g. momentum or a chosen enthalpy datum.
    pub carried: [[f64; N]; 2],
}

impl<const N: usize> Default for EdgeDensityTransfer<N> {
    fn default() -> Self {
        Self {
            mass_kg: [0.0; 2],
            carried: [[0.0; N]; 2],
        }
    }
}

/// Physical inventories and roundoff balance defects of one accepted step.
/// Enthalpy transport alone does not include pressure-work or thermodynamic
/// closure; momentum mixing loss is not deposited as heat by this operator.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityTransportReport<const N: usize> {
    /// Accepted interval [s].
    pub timestep_s: f64,
    /// Actual largest outgoing volume fraction.
    pub maximum_courant: f64,
    /// Domain mass before and after [kg].
    pub mass_before_kg: f64,
    pub mass_after_kg: f64,
    /// Integrated physical boundary influx and efflux [kg].
    pub inlet_mass_kg: f64,
    pub outlet_mass_kg: f64,
    /// Sum of cell mass increments plus efflux minus influx [kg].
    pub mass_balance_residual_kg: f64,
    /// Mass residual divided by the largest inventory or boundary mass scale.
    pub relative_mass_balance_residual: f64,
    /// Extensive domain inventories before and after.
    pub carried_before: [f64; N],
    pub carried_after: [f64; N],
    /// Extensive boundary influx and efflux in the declared component units.
    pub carried_in: [f64; N],
    pub carried_out: [f64; N],
    /// Direct cell increments plus efflux minus influx; not a conservation repair.
    pub carried_balance_residual: [f64; N],
    /// Signed residual / maximum absolute inventory or boundary scale.
    pub relative_carried_balance_residual: [f64; N],
}

/// Invalid geometry/state, excessive work storage/CFL, or interrupted transport.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DensityTransportError {
    /// The original BDM face evaluation refused its input.
    Volume(FluxTransportError),
    /// Dimensions, values or topology do not match the declared model.
    Invalid(&'static str),
    /// Prescribed volume flow is not locally divergence-free within allowance.
    VolumeDivergence { cell: usize, relative_residual: f64 },
    /// The requested interval violates the admitted outgoing-volume CFL limit.
    Courant { actual: f64, maximum: f64 },
    /// Required logical vector storage exceeds its allowance.
    Storage { required: usize, maximum: usize },
    /// A fallible allocation could not reserve its requested storage.
    Allocation,
    /// A required physical inventory, rate, or update is not representable.
    Arithmetic(&'static str),
    /// No output or accepted face-transfer view was published.
    Cancelled,
}

impl core::fmt::Display for DensityTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "BDM density transport: {self:?}")
    }
}
impl std::error::Error for DensityTransportError {}
impl From<FluxTransportError> for DensityTransportError {
    fn from(value: FluxTransportError) -> Self {
        Self::Volume(value)
    }
}
type Result<T> = core::result::Result<T, DensityTransportError>;

/// Original fixed geometry with a frozen, admitted BDM volume velocity.
/// Reprepare when that velocity changes. No material source is queried here.
pub struct PreparedDensityTransport<'m> {
    pub(super) mesh: &'m TriMesh,
    volumes: Vec<f64>,
    volume_rates: Vec<EdgeVolumeTransport>,
    // Retain the original affine normal traces for compatible momentum flux.
    pub(super) velocity_moments: Vec<[f64; 2]>,
    // Per-direction rate divided by its donor area; a boundary inlet uses
    // the receiver area to avoid constructing a fictitious exterior cell.
    donor_rates: Vec<[f64; 2]>,
    maximum_outgoing_rate_s: f64,
    maximum_relative_volume_divergence: f64,
}

/// Caller-owned reusable scratch. A warm step allocates no heap storage.
/// Its contents may change on refusal, but `face_transfers` returns None until
/// an entire step succeeds. Accepted output cells are copied only at the end.
pub struct DensityTransportWorkspace<const N: usize> {
    candidate: Vec<DensityTransportCell<N>>,
    masses: Vec<f64>,
    mass_change: Vec<Sum>,
    specific_change: Vec<[Sum; N]>,
    transfers: Vec<EdgeDensityTransfer<N>>,
    accepted: bool,
}

impl<const N: usize> DensityTransportWorkspace<N> {
    /// Both directions and their actual carried amounts from the latest
    /// successful step; invalidated at the start of every subsequent attempt.
    pub fn face_transfers(&self) -> Option<&[EdgeDensityTransfer<N>]> {
        self.accepted.then_some(self.transfers.as_slice())
    }
}

impl<'m> FluxSystem<'m> {
    /// Bind conservative density transport to this original mesh and velocity.
    /// `depth_m` converts 2-D areas/rates to physical volumes/mass. The per-cell
    /// relative volume-divergence allowance lies in [0,1e-12] and is numerical
    /// only; no prescribed velocity or face rate is projected or repaired.
    /// The byte cap includes constructor scratch plus retained vector elements.
    /// Cancellation is checked per original cell/edge, before any publication.
    pub fn prepare_density_transport(
        &self,
        velocity: &[f64],
        depth_m: f64,
        relative_volume_divergence_tolerance: f64,
        max_storage_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedDensityTransport<'m>> {
        poll(&mut cancelled)?;
        let mesh = self.mesh;
        let nc = mesh.tris.len();
        let ne = mesh.edges.len();
        if nc == 0
            || mesh.areas.len() != nc
            || mesh.tri_edges.len() != nc
            || ne.checked_mul(2) != Some(self.n_u)
            || !positive(depth_m)
            || !relative_volume_divergence_tolerance.is_finite()
            || !(0.0..=1e-12).contains(&relative_volume_divergence_tolerance)
        {
            return Err(DensityTransportError::Invalid(
                "positive depth, original nonempty mesh and bounded divergence tolerance required",
            ));
        }
        if velocity.len() != self.n_u {
            return Err(FluxTransportError::VelocityLength {
                expected: self.n_u,
                actual: velocity.len(),
            }
            .into());
        }
        admit_storage(
            &[
                (nc, size_of::<f64>() + 2 * size_of::<Sum>()),
                (ne, size_of::<EdgeVolumeTransport>() + 2 * size_of::<[f64; 2]>()),
            ],
            max_storage_bytes,
        )?;
        let mut volumes = allocate(nc, 0.0)?;
        let mut incoming = allocate(nc, Sum::default())?;
        let mut outgoing = allocate(nc, Sum::default())?;
        let mut volume_rates = allocate(ne, EdgeVolumeTransport::default())?;
        let mut velocity_moments = allocate(ne, [0.0; 2])?;
        let mut donor_rates = allocate(ne, [0.0; 2])?;
        for (cell, &area) in mesh.areas.iter().enumerate() {
            poll(&mut cancelled)?;
            if !positive(area) {
                return Err(DensityTransportError::Invalid(
                    "finite positive original cell areas required",
                ));
            }
            volumes[cell] = positive_value(area * depth_m, "cell volume")?;
            let faces = mesh.tri_edges[cell];
            for (local, &(edge, sign)) in faces.iter().enumerate() {
                let adjacent = mesh
                    .edges
                    .get(edge)
                    .ok_or(DensityTransportError::Invalid("cell edge index"))?
                    .tris;
                if faces[..local].iter().any(|&(previous, _)| previous == edge)
                    || !((adjacent.0 == cell && sign == 1.0)
                        || (adjacent.1 == cell && sign == -1.0))
                {
                    return Err(DensityTransportError::Invalid(
                        "original oriented cell edges required",
                    ));
                }
            }
        }
        for (index, edge) in mesh.edges.iter().enumerate() {
            poll(&mut cancelled)?;
            let (i, j) = edge.tris;
            if i >= nc || (j != usize::MAX && (j >= nc || i == j)) {
                return Err(DensityTransportError::Invalid(
                    "original edge adjacency required",
                ));
            }
            if !mesh.tri_edges[i].contains(&(index, 1.0))
                || (j != usize::MAX && !mesh.tri_edges[j].contains(&(index, -1.0)))
            {
                return Err(DensityTransportError::Invalid(
                    "every global edge must belong to its original oriented cell faces",
                ));
            }
            for dof in [2 * index, 2 * index + 1] {
                if !velocity[dof].is_finite() {
                    return Err(FluxTransportError::NonFiniteVelocity { dof }.into());
                }
            }
            let rate = edge_volume_rate(
                index,
                edge.len,
                [velocity[2 * index], velocity[2 * index + 1]],
            )?;
            let forward = rate.along_normal_m2_s();
            let reverse = rate.against_normal_m2_s();
            outgoing[i].add(forward)?;
            incoming[i].add(reverse)?;
            let reverse_cell = if j == usize::MAX {
                i
            } else {
                outgoing[j].add(reverse)?;
                incoming[j].add(forward)?;
                j
            };
            donor_rates[index] = [
                quotient(forward, mesh.areas[i])?,
                quotient(reverse, mesh.areas[reverse_cell])?,
            ];
            volume_rates[index] = rate;
            velocity_moments[index] = [velocity[2 * index], velocity[2 * index + 1]];
        }
        let mut maximum_outgoing_rate_s: f64 = 0.0;
        let mut maximum_relative_volume_divergence: f64 = 0.0;
        for cell in 0..nc {
            poll(&mut cancelled)?;
            let incoming = incoming[cell].value()?;
            let outgoing = outgoing[cell].value()?;
            let scale = incoming.max(outgoing);
            let relative = if scale == 0.0 {
                0.0
            } else {
                (incoming - outgoing).abs() / scale
            };
            if relative > relative_volume_divergence_tolerance {
                return Err(DensityTransportError::VolumeDivergence {
                    cell,
                    relative_residual: relative,
                });
            }
            maximum_relative_volume_divergence = maximum_relative_volume_divergence.max(relative);
            maximum_outgoing_rate_s =
                maximum_outgoing_rate_s.max(quotient(outgoing, mesh.areas[cell])?);
        }
        poll(&mut cancelled)?;
        Ok(PreparedDensityTransport {
            mesh,
            volumes,
            volume_rates,
            velocity_moments,
            donor_rates,
            maximum_outgoing_rate_s,
            maximum_relative_volume_divergence,
        })
    }
}

impl PreparedDensityTransport<'_> {
    /// Physical volumes [m³] in original triangle order.
    pub fn cell_volumes_m3(&self) -> &[f64] {
        &self.volumes
    }
    /// The same native sign-split volume rates, in original edge order [m²/s].
    pub fn edge_volume_rates(&self) -> &[EdgeVolumeTransport] {
        &self.volume_rates
    }
    /// Actual divergence left in the prescribed field, not a repaired value.
    pub const fn maximum_relative_volume_divergence(&self) -> f64 {
        self.maximum_relative_volume_divergence
    }

    /// Largest admissible interval [s] for a Courant cap strictly in (0,1).
    /// Returns infinity at rest or when no finite positive timestep exceeds it.
    pub fn maximum_timestep_s(&self, maximum_courant: f64) -> Result<f64> {
        admit_courant(maximum_courant)?;
        Ok(maximum_courant / self.maximum_outgoing_rate_s)
    }

    /// Allocate bounded scratch once for N specific quantities. This cap covers
    /// scratch vector elements, excluding the plan, inputs and output buffers.
    pub fn workspace<const N: usize>(
        &self,
        max_storage_bytes: usize,
    ) -> Result<DensityTransportWorkspace<N>> {
        let nc = self.volumes.len();
        let ne = self.volume_rates.len();
        admit_storage(
            &[
                (
                    nc,
                    size_of::<DensityTransportCell<N>>()
                        + size_of::<f64>()
                        + size_of::<Sum>()
                        + size_of::<[Sum; N]>(),
                ),
                (ne, size_of::<EdgeDensityTransfer<N>>()),
            ],
            max_storage_bytes,
        )?;
        Ok(DensityTransportWorkspace {
            candidate: allocate(nc, DensityTransportCell::default())?,
            masses: allocate(nc, 0.0)?,
            mass_change: allocate(nc, Sum::default())?,
            specific_change: allocate(nc, [Sum::default(); N])?,
            transfers: allocate(ne, EdgeDensityTransfer::default())?,
            accepted: false,
        })
    }

    /// Advance density and all carried inventories synchronously from `cells`.
    /// `inflow` has one entry per original edge, with an explicit positive state
    /// on every boundary edge carrying reverse flow; interior entries are None.
    /// Outgoing boundary transport uses the original interior donor, including
    /// on an edge which ALSO carries an incoming stream.
    ///
    /// With m=rho*V and each incoming transfer a, the updates are
    /// `m_new=m_old+sum(in)-sum(out)` and
    /// `s_new=s_old+sum_in(a*(s_donor-s_old))/m_new`.
    /// They conserve m and m*s up to reported floating-point roundoff, preserve
    /// constant specific quantities exactly, and preserve positivity under the
    /// outgoing-volume CFL limit. Density bounds additionally use the admitted
    /// local volume balance; numerical divergence is not erased.
    ///
    /// Work is O((cells+edges)*N), with no warm heap allocation. Scratch can be
    /// reused across equal-sized plans. Output is unchanged on any refusal or
    /// cancellation; scratch's accepted-transfer view is then invalidated.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_into<const N: usize>(
        &self,
        cells: &[DensityTransportCell<N>],
        inflow: &[Option<DensityTransportCell<N>>],
        timestep_s: f64,
        maximum_courant: f64,
        output: &mut [DensityTransportCell<N>],
        scratch: &mut DensityTransportWorkspace<N>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DensityTransportReport<N>> {
        scratch.accepted = false;
        poll(&mut cancelled)?;
        let nc = self.volumes.len();
        let ne = self.volume_rates.len();
        if cells.len() != nc
            || output.len() != nc
            || inflow.len() != ne
            || scratch.candidate.len() != nc
            || scratch.transfers.len() != ne
            || !positive(timestep_s)
        {
            return Err(DensityTransportError::Invalid(
                "matching cell/edge/workspace dimensions and positive finite timestep required",
            ));
        }
        admit_courant(maximum_courant)?;
        let courant = finite(timestep_s * self.maximum_outgoing_rate_s, "Courant number")?;
        if courant > maximum_courant {
            return Err(DensityTransportError::Courant {
                actual: courant,
                maximum: maximum_courant,
            });
        }
        let mut mass_before = Sum::default();
        let mut carried_before = [Sum::default(); N];
        let mut absolute_before = [Sum::default(); N];
        for (cell, state) in cells.iter().enumerate() {
            poll(&mut cancelled)?;
            admit_cell(state)?;
            let mass = positive_value(state.density_kg_m3 * self.volumes[cell], "cell mass")?;
            scratch.masses[cell] = mass;
            scratch.mass_change[cell] = Sum::default();
            scratch.specific_change[cell] = [Sum::default(); N];
            mass_before.add(mass)?;
            for component in 0..N {
                let value = finite(
                    mass * state.specific[component],
                    "initial carried inventory",
                )?;
                carried_before[component].add(value)?;
                absolute_before[component].add(value.abs())?;
            }
        }
        let mut inlet_mass = Sum::default();
        let mut outlet_mass = Sum::default();
        let mut carried_in = [Sum::default(); N];
        let mut carried_out = [Sum::default(); N];
        let mut absolute_boundary = [Sum::default(); N];
        let mut mass_residual = Sum::default();
        let mut carried_residual = [Sum::default(); N];
        for (index, edge) in self.mesh.edges.iter().enumerate() {
            poll(&mut cancelled)?;
            let (i, j) = edge.tris;
            if j != usize::MAX && inflow[index].is_some() {
                return Err(DensityTransportError::Invalid(
                    "inflow states belong only to boundary edges",
                ));
            }
            if let Some(state) = &inflow[index] {
                admit_cell(state)?;
            }
            let mut transfer = EdgeDensityTransfer::default();
            for direction in 0..2 {
                let rate = self.donor_rates[index][direction];
                if rate == 0.0 {
                    continue;
                }
                let fraction = positive_value(timestep_s * rate, "directional volume fraction")?;
                let (donor, receiver) = if direction == 0 { (i, j) } else { (j, i) };
                let state = if donor == usize::MAX {
                    inflow[index]
                        .as_ref()
                        .ok_or(DensityTransportError::Invalid(
                            "every incoming boundary stream needs a physical donor state",
                        ))?
                } else {
                    &cells[donor]
                };
                let mass = if donor == usize::MAX {
                    dimensional_rate(fraction, state.density_kg_m3, self.volumes[receiver])
                        .ok_or(DensityTransportError::Arithmetic("inlet mass transfer"))?
                } else {
                    positive_value(fraction * scratch.masses[donor], "donor mass transfer")?
                };
                transfer.mass_kg[direction] = mass;
                if donor != usize::MAX {
                    scratch.mass_change[donor].add(-mass)?;
                } else {
                    inlet_mass.add(mass)?;
                    mass_residual.add(-mass)?;
                }
                if receiver != usize::MAX {
                    scratch.mass_change[receiver].add(mass)?;
                } else {
                    outlet_mass.add(mass)?;
                    mass_residual.add(mass)?;
                }
                for component in 0..N {
                    let value = finite(mass * state.specific[component], "carried face transfer")?;
                    transfer.carried[direction][component] = value;
                    if receiver != usize::MAX {
                        let local = cells[receiver].specific[component];
                        let difference = state.specific[component] - local;
                        // The weighted difference can be representable even if
                        // opposite finite specific values overflow subtraction.
                        let change = if difference.is_finite() {
                            mass * difference
                        } else {
                            value - mass * local
                        };
                        scratch.specific_change[receiver][component]
                            .add(finite(change, "incoming specific-inventory increment")?)?;
                    }
                    if donor == usize::MAX {
                        carried_in[component].add(value)?;
                        carried_residual[component].add(-value)?;
                        absolute_boundary[component].add(value.abs())?;
                    } else if receiver == usize::MAX {
                        carried_out[component].add(value)?;
                        carried_residual[component].add(value)?;
                        absolute_boundary[component].add(value.abs())?;
                    }
                }
            }
            scratch.transfers[index] = transfer;
        }
        let mut mass_after = Sum::default();
        let mut carried_after = [Sum::default(); N];
        let mut absolute_after = [Sum::default(); N];
        for (cell, state) in cells.iter().enumerate() {
            poll(&mut cancelled)?;
            let old_mass = scratch.masses[cell];
            let new_mass = positive_value(
                old_mass + scratch.mass_change[cell].value()?,
                "positive updated cell mass",
            )?;
            let density = positive_value(new_mass / self.volumes[cell], "updated density")?;
            // Report the inventory represented by the output density, including
            // its conversion roundoff rather than a hidden higher-precision mass.
            let represented_mass =
                positive_value(density * self.volumes[cell], "represented updated mass")?;
            let mut specific = state.specific;
            mass_after.add(represented_mass)?;
            mass_residual.add(represented_mass - old_mass)?;
            for component in 0..N {
                let change = scratch.specific_change[cell][component].value()?;
                if change != 0.0 {
                    let candidate = specific[component] + change / new_mass;
                    specific[component] = if candidate.is_finite()
                        && candidate.abs() >= 0.5 * specific[component].abs()
                    {
                        candidate
                    } else {
                        // A dense incoming stream with nearly zero specific
                        // value must not cancel away the retained donor's
                        // finite inventory. The original three faces provide
                        // the direct convex mixture without additional storage.
                        self.extensive_mix(
                            cell,
                            component,
                            old_mass,
                            state.specific[component],
                            new_mass,
                            &scratch.transfers,
                        )?
                    };
                }
                let value = finite(
                    represented_mass * specific[component],
                    "updated carried inventory",
                )?;
                carried_after[component].add(value)?;
                absolute_after[component].add(value.abs())?;
                // Compensated product differences audit the represented state
                // without cancellation of either a large common datum or a
                // large density/specific contrast. No physical state is repaired.
                let old_value = old_mass * state.specific[component];
                carried_residual[component].add(value)?;
                carried_residual[component].add(-old_value)?;
                carried_residual[component]
                    .add(represented_mass.mul_add(specific[component], -value))?;
                carried_residual[component]
                    .add(-old_mass.mul_add(state.specific[component], -old_value))?;
            }
            scratch.candidate[cell] = DensityTransportCell {
                density_kg_m3: density,
                specific,
            };
        }
        let mass_before = mass_before.value()?;
        let mass_after = mass_after.value()?;
        let inlet_mass = inlet_mass.value()?;
        let outlet_mass = outlet_mass.value()?;
        let mass_residual = mass_residual.value()?;
        let mut report = DensityTransportReport {
            timestep_s,
            maximum_courant: courant,
            mass_before_kg: mass_before,
            mass_after_kg: mass_after,
            inlet_mass_kg: inlet_mass,
            outlet_mass_kg: outlet_mass,
            mass_balance_residual_kg: mass_residual,
            relative_mass_balance_residual: mass_residual
                / mass_before.max(mass_after).max(inlet_mass).max(outlet_mass),
            carried_before: [0.0; N],
            carried_after: [0.0; N],
            carried_in: [0.0; N],
            carried_out: [0.0; N],
            carried_balance_residual: [0.0; N],
            relative_carried_balance_residual: [0.0; N],
        };
        for component in 0..N {
            report.carried_before[component] = carried_before[component].value()?;
            report.carried_after[component] = carried_after[component].value()?;
            report.carried_in[component] = carried_in[component].value()?;
            report.carried_out[component] = carried_out[component].value()?;
            let residual = carried_residual[component].value()?;
            report.carried_balance_residual[component] = residual;
            let scale = absolute_before[component]
                .value()?
                .max(absolute_after[component].value()?)
                .max(absolute_boundary[component].value()?);
            report.relative_carried_balance_residual[component] =
                if scale == 0.0 { 0.0 } else { residual / scale };
        }
        poll(&mut cancelled)?;
        output.copy_from_slice(&scratch.candidate);
        scratch.accepted = true;
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    fn extensive_mix<const N: usize>(
        &self,
        cell: usize,
        component: usize,
        old_mass: f64,
        specific: f64,
        new_mass: f64,
        transfers: &[EdgeDensityTransfer<N>],
    ) -> Result<f64> {
        let mut outgoing = Sum::default();
        let mut inventory = Sum::default();
        for &(edge, sign) in &self.mesh.tri_edges[cell] {
            let out = usize::from(sign < 0.0);
            outgoing.add(transfers[edge].mass_kg[out])?;
            inventory.add(transfers[edge].carried[1 - out][component])?;
        }
        let retained = positive_value(old_mass - outgoing.value()?, "retained donor mass")?;
        inventory.add(finite(retained * specific, "retained carried inventory")?)?;
        finite(inventory.value()? / new_mass, "mixed specific quantity")
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<()> {
        let next = finite(self.value + value, "inventory reduction")?;
        self.correction = finite(
            self.correction
                + if self.value.abs() >= value.abs() {
                    (self.value - next) + value
                } else {
                    (value - next) + self.value
                },
            "inventory compensation",
        )?;
        self.value = next;
        Ok(())
    }
    fn value(self) -> Result<f64> {
        finite(self.value + self.correction, "inventory sum")
    }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: f64, context: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DensityTransportError::Arithmetic(context))
    }
}
fn positive_value(value: f64, context: &'static str) -> Result<f64> {
    if positive(value) {
        Ok(value)
    } else {
        Err(DensityTransportError::Arithmetic(context))
    }
}
fn quotient(rate: f64, area: f64) -> Result<f64> {
    if rate == 0.0 {
        Ok(0.0)
    } else {
        positive_value(rate / area, "volume rate divided by cell area")
    }
}
fn admit_cell<const N: usize>(state: &DensityTransportCell<N>) -> Result<()> {
    if positive(state.density_kg_m3) && state.specific.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(DensityTransportError::Invalid(
            "positive density and finite specific quantities required",
        ))
    }
}
fn admit_courant(value: f64) -> Result<()> {
    if value.is_finite() && value > 0.0 && value < 1.0 {
        Ok(())
    } else {
        Err(DensityTransportError::Invalid(
            "maximum Courant number must lie strictly in (0,1)",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(DensityTransportError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit_storage(parts: &[(usize, usize)], maximum: usize) -> Result<()> {
    let required = parts
        .iter()
        .try_fold(0usize, |sum, &(count, size)| {
            sum.checked_add(count.checked_mul(size)?)
        })
        .ok_or(DensityTransportError::Arithmetic("storage size"))?;
    if required > maximum {
        Err(DensityTransportError::Storage { required, maximum })
    } else {
        Ok(())
    }
}
fn allocate<T: Clone>(count: usize, value: T) -> Result<Vec<T>> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| DensityTransportError::Allocation)?;
    result.resize(count, value);
    Ok(result)
}
