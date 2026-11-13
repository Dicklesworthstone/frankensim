//! One original triangle inventory and exact directional BDM transport map.

use std::collections::BTreeMap;

use fs_flux::{
    FluxSystem, MaterialNavierStokesBoundary, MaterialStokesBoundary, MaterialStokesSolution,
    StokesPressureReference,
};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::advection::{Advected, MassBoundary, MassFlow};
use fs_thermal::enthalpy::backward_euler::{
    BoundaryStepCoefficients, StepCoefficients, TransportStepCoefficients,
};
use fs_thermal::enthalpy::{CellGeometry, ThermalStorageBasis};

use super::*;
use crate::stokes_motion::heating::StokesHeatBinding;

struct Storage {
    material: MaterialId,
    mass_kg: f64,
    geometry: CellGeometry,
    caloric_identity: [u8; 32],
}

pub(super) struct ThermalBinding {
    cells: Vec<CellId>,
    storage: Vec<Storage>,
    plane_z_m: f64,
    depth_m: f64,
    density_kg_m3: f64,
    maximum_connections: usize,
}

impl ThermalBinding {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        cx: &Cx<'_>,
        system: &FluxSystem<'_>,
        thermal: &DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        depth_m: f64,
        density_kg_m3: f64,
        maximum_connections: usize,
    ) -> Result<Self, LaminarEnthalpyError> {
        if !depth_m.is_finite()
            || depth_m <= 0.0
            || !density_kg_m3.is_finite()
            || density_kg_m3 <= 0.0
            || cells.triangle_cells.is_empty()
            || system
                .mesh
                .edges
                .len()
                .checked_mul(2)
                .is_none_or(|n| n > maximum_connections)
        {
            return Err(LaminarEnthalpyError::Invalid(
                "positive extrusion/density and a complete bounded two-direction face map required",
            ));
        }
        StokesHeatBinding {
            mesh: system.mesh,
            density_kg_m3,
            depth_m,
            time_s: thermal.time_s(),
        }
        .check_thermal_cells(thermal, cells, &mut || cx.checkpoint().is_err())?;
        let mut storage = Vec::with_capacity(cells.triangle_cells.len());
        let material = thermal.cells()[&cells.triangle_cells[0]].material_id;
        for &id in cells.triangle_cells {
            checkpoint(cx)?;
            let cell = &thermal.cells()[&id];
            let caloric_identity = thermal.material_caloric_identity(cell.material_id)?.ok_or(
                LaminarEnthalpyError::Invalid(
                    "an original identity-bearing fluid caloric law is required",
                ),
            )?;
            if cell.material_id != material
                || thermal.material_storage_basis(cell.material_id)?
                    != ThermalStorageBasis::Enthalpy
            {
                return Err(LaminarEnthalpyError::Invalid(
                    "all fluid cells must share the same original enthalpy material",
                ));
            }
            storage.push(Storage {
                material,
                mass_kg: cell.mass_kg,
                geometry: cell.geometry,
                caloric_identity,
            });
        }
        Ok(Self {
            cells: cells.triangle_cells.to_vec(),
            storage,
            plane_z_m: cells.plane_z_m,
            depth_m,
            density_kg_m3,
            maximum_connections,
        })
    }

    pub(super) fn cells(&self) -> &[CellId] {
        &self.cells
    }

    pub(super) const fn depth_m(&self) -> f64 {
        self.depth_m
    }

    fn native<'a>(&self, system: &'a FluxSystem<'_>, time_s: f64) -> StokesHeatBinding<'a> {
        StokesHeatBinding {
            mesh: system.mesh,
            density_kg_m3: self.density_kg_m3,
            depth_m: self.depth_m,
            time_s,
        }
    }

    pub(super) fn check<E>(
        &self,
        cx: &Cx<'_>,
        system: &FluxSystem<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), LaminarEnthalpyError<E>> {
        self.native(system, thermal.time_s()).check_thermal_cells(
            thermal,
            StokesThermalCells {
                triangle_cells: &self.cells,
                plane_z_m: self.plane_z_m,
            },
            &mut || cx.checkpoint().is_err(),
        )?;
        for (&id, original) in self.cells.iter().zip(&self.storage) {
            checkpoint(cx)?;
            let cell = &thermal.cells()[&id];
            if cell.material_id != original.material
                || cell.mass_kg != original.mass_kg
                || cell.geometry != original.geometry
                || thermal.material_caloric_identity(cell.material_id)?
                    != Some(original.caloric_identity)
                || thermal.material_storage_basis(cell.material_id)?
                    != ThermalStorageBasis::Enthalpy
            {
                return Err(LaminarEnthalpyError::Invalid(
                    "original fluid inventory, geometry or caloric identity changed",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn temperatures<E>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<Vec<f64>, LaminarEnthalpyError<E>> {
        self.cells
            .iter()
            .map(|&id| {
                checkpoint(cx)?;
                thermal.temperature_k(id).map_err(Into::into)
            })
            .collect()
    }

    pub(super) fn slip_routes<E>(
        &self,
        cx: &Cx<'_>,
        system: &FluxSystem<'_>,
        thermal: &DepositionSimulation<'_>,
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        routes: &[StokesSlipHeatRoute],
    ) -> Result<BTreeMap<usize, StokesSlipHeatRoute>, LaminarEnthalpyError<E>> {
        let native: Vec<_> = boundaries
            .iter()
            .map(|boundary| match *boundary {
                MaterialNavierStokesBoundary::Velocity { edge, value } => {
                    MaterialStokesBoundary::Velocity { edge, value }
                }
                MaterialNavierStokesBoundary::NavierSlip {
                    edge,
                    value,
                    tangential_resistance_pa_s_per_m,
                } => MaterialStokesBoundary::NavierSlip {
                    edge,
                    value,
                    tangential_resistance_pa_s_per_m,
                },
                MaterialNavierStokesBoundary::Open { edge, traction, .. }
                | MaterialNavierStokesBoundary::OpenReservoir { edge, traction, .. } => {
                    MaterialStokesBoundary::TotalTraction {
                        edge,
                        value: traction,
                    }
                }
            })
            .collect();
        self.native(system, thermal.time_s())
            .check_slip_heat_routes(
                thermal,
                StokesThermalCells {
                    triangle_cells: &self.cells,
                    plane_z_m: self.plane_z_m,
                },
                &[],
                &native,
                routes,
                &mut || cx.checkpoint().is_err(),
            )
            .map_err(Into::into)
    }

    // The original signed normal trace can carry two different donors on the
    // same face. Never replace these rates with an absolute/net mean flux.
    pub(super) fn advection<E>(
        &self,
        cx: &Cx<'_>,
        system: &FluxSystem<'_>,
        velocity: &[f64],
        incoming: &[Option<LaminarInletEnthalpy>],
    ) -> Result<Advected<()>, LaminarEnthalpyError<E>> {
        if incoming.len() != system.mesh.edges.len() {
            return Err(LaminarEnthalpyError::Invalid(
                "one optional inlet state per original edge required",
            ));
        }
        let rates = system
            .edge_volume_transport(velocity)
            .map_err(LaminarEnthalpyError::FaceTransport)?;
        let mut flows = Vec::new();
        let mut boundaries = Vec::new();
        for (edge, (geometry, rate)) in system.mesh.edges.iter().zip(&rates).enumerate() {
            checkpoint(cx)?;
            let id = u64::try_from(edge)
                .ok()
                .and_then(|n| n.checked_mul(2))
                .ok_or(LaminarEnthalpyError::Invalid(
                    "native face identity overflow",
                ))?;
            let outgoing = product(&[rate.along_normal_m2_s(), self.depth_m, self.density_kg_m3])?;
            let ingoing = product(&[rate.against_normal_m2_s(), self.depth_m, self.density_kg_m3])?;
            let a = self.cells[geometry.tris.0];
            if system.mesh.is_boundary(edge) {
                if incoming[edge].is_some_and(|input| !input.specific_h_j_kg.is_finite()) {
                    return Err(LaminarEnthalpyError::Invalid(
                        "finite incoming specific enthalpy required",
                    ));
                }
                if outgoing > 0.0 {
                    boundaries.push(MassBoundary::Outlet {
                        id,
                        cell: a,
                        kg_s: outgoing,
                    });
                }
                if ingoing > 0.0 {
                    let input = incoming[edge].ok_or(LaminarEnthalpyError::Invalid(
                        "every actually incoming face needs its original external enthalpy",
                    ))?;
                    boundaries.push(MassBoundary::Inlet {
                        id: id + 1,
                        cell: a,
                        kg_s: ingoing,
                        material_id: input.material_id,
                        specific_h_j_kg: input.specific_h_j_kg,
                    });
                }
            } else {
                if incoming[edge].is_some() {
                    return Err(LaminarEnthalpyError::Invalid(
                        "interior faces cannot carry external inlet states",
                    ));
                }
                let b = self.cells[geometry.tris.1];
                if outgoing > 0.0 {
                    flows.push(MassFlow {
                        id,
                        a,
                        b,
                        kg_s: outgoing,
                    });
                }
                if ingoing > 0.0 {
                    flows.push(MassFlow {
                        id: id + 1,
                        a: b,
                        b: a,
                        kg_s: ingoing,
                    });
                }
            }
        }
        Ok(Advected::new((), flows, self.maximum_connections)?
            .with_boundaries(boundaries, self.maximum_connections)?)
    }

    pub(super) fn coefficients<T, E>(
        &self,
        cx: &Cx<'_>,
        system: &FluxSystem<'_>,
        thermal: &DepositionSimulation<'_>,
        flow: &MaterialStokesSolution,
        routes: &BTreeMap<usize, StokesSlipHeatRoute>,
        response: LaminarHeatExchange<T>,
    ) -> Result<(TransportStepCoefficients<T>, [f64; 4]), LaminarEnthalpyError<E>> {
        let mut advected = self
            .advection::<E>(cx, system, flow.velocity_dofs(), &response.inlet_enthalpies)?
            .endpoint_transport(thermal, || thermal_checkpoint(cx))?;
        let other_power = thermal.external_step_power_with_boundary_conductances(
            &response.heat.bulk.sources_w,
            &response.heat.bulk.conductances_w_k,
            &response.heat.boundaries,
            || thermal_checkpoint(cx),
        )?;
        let external_stream = checked_sum(
            &advected
                .heat
                .bulk
                .sources_w
                .values()
                .copied()
                .collect::<Vec<_>>(),
        )?;
        let BoundaryStepCoefficients { bulk, boundaries } = response.heat;
        let StepCoefficients {
            mut sources_w,
            conductances_w_k,
            payload,
        } = bulk;
        let heating = flow.viscous_dissipation().cell_heating_w_m3();
        if heating.len() != self.cells.len() {
            return Err(LaminarEnthalpyError::Invalid(
                "incomplete accepted native strain heat",
            ));
        }
        let mut bulk_power = Dd::from_f64(0.0);
        for (cell, (&id, &rate)) in self.cells.iter().zip(heating).enumerate() {
            checkpoint(cx)?;
            let power = product(&[rate, system.mesh.areas[cell], self.depth_m])?;
            add_source(&mut sources_w, id, power)?;
            bulk_power = bulk_power + Dd::from_f64(power);
        }
        let mut slip_power = Dd::from_f64(0.0);
        if !routes.is_empty() {
            // Slip heat is pressure-gauge independent. Querying the existing
            // reaction owner uses zero for an arbitrary mean-zero gauge, and
            // preserves every traction datum. No force/absolute pressure is
            // inferred or returned from this auxiliary loss query.
            let mut gauges = vec![None; system.pressure_component_count()];
            for (cell, &component) in system.pressure_components().iter().enumerate() {
                gauges[component] = match flow.pressure_reference_at(cell) {
                    Some(StokesPressureReference::MeanZero) => Some(0.0),
                    Some(StokesPressureReference::TractionDatum) => None,
                    _ => {
                        return Err(LaminarEnthalpyError::Invalid(
                            "missing local native pressure reference",
                        ));
                    }
                };
            }
            for (&edge, route) in routes {
                checkpoint(cx)?;
                let loss = flow
                    .boundary_loads(&[edge], [0.0; 2], &gauges)
                    .map_err(LaminarEnthalpyError::Momentum)?
                    .slip_dissipation_w_per_m;
                let power = product(&[loss, self.depth_m])?;
                let fraction = route.fluid_heat_fraction;
                let (fluid, wall) = if fraction <= 0.5 {
                    let fluid = product(&[power, fraction])?;
                    (fluid, power - fluid)
                } else {
                    let wall = product(&[power, 1.0 - fraction])?;
                    (power - wall, wall)
                };
                add_source(
                    &mut sources_w,
                    self.cells[system.mesh.edges[edge].tris.0],
                    fluid,
                )?;
                add_source(&mut sources_w, route.wall_cell, wall)?;
                slip_power = slip_power + Dd::from_f64(power);
            }
        }
        for (id, watts) in std::mem::take(&mut advected.heat.bulk.sources_w) {
            add_source(&mut sources_w, id, watts)?;
        }
        let powers = [
            bulk_power.to_f64(),
            slip_power.to_f64(),
            external_stream,
            other_power,
        ];
        if powers.iter().any(|value| !value.is_finite()) {
            return Err(DepositionError::Arithmetic.into());
        }
        Ok((
            TransportStepCoefficients {
                heat: BoundaryStepCoefficients {
                    bulk: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload,
                    },
                    boundaries,
                },
                transfers: advected.transfers,
            },
            powers,
        ))
    }
}

pub(super) fn product<const N: usize, E>(
    factors: &[f64; N],
) -> Result<f64, LaminarEnthalpyError<E>> {
    if factors
        .iter()
        .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(DepositionError::Arithmetic.into());
    }
    if factors.contains(&0.0) {
        return Ok(0.0);
    }
    let mut sorted = *factors;
    sorted.sort_unstable_by(f64::total_cmp);
    let mut value = 1.0;
    let (mut low, mut high) = (0, N);
    while high - low >= 2 {
        high -= 1;
        value *= sorted[low] * sorted[high];
        low += 1;
    }
    if low < high {
        value *= sorted[low];
    }
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}

fn add_source<E>(
    sources: &mut BTreeMap<CellId, f64>,
    id: CellId,
    watts: f64,
) -> Result<(), LaminarEnthalpyError<E>> {
    let total = sources.entry(id).or_insert(0.0);
    *total = checked_sum(&[*total, watts])?;
    Ok(())
}

pub(super) fn checked_sum<E>(terms: &[f64]) -> Result<f64, LaminarEnthalpyError<E>> {
    let mut sum = Dd::from_f64(0.0);
    for &value in terms {
        if !value.is_finite() {
            return Err(DepositionError::Arithmetic.into());
        }
        sum = sum + Dd::from_f64(value);
    }
    let result = sum.to_f64();
    if result.is_finite() {
        Ok(result)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
