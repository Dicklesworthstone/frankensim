//! One source constitutive root for molecular plus modeled shear viscosity.
//!
//! Smagorinsky depends on the trial deviatoric rate; WALE/Vreman use the actual
//! full liquid gradient. Both close gamma*tau_total = Guo stress numerator.
//! The same material response supplies independent bulk viscosity. Native
//! collision, mass transport and conversion remain the only evolution owner.

use std::borrow::Cow;

use super::*;
use crate::d3q19::freesurface3::FreeSurfaceKinematicsWorkspace3;
use crate::d3q19::freesurface3::sgs::{
    FreeSurfaceSgsError3, FreeSurfaceSgsReport3, FreeSurfaceSgsWorkspace3, gradient_sgs_at,
    molecular_relaxation,
};
use crate::lattice3::{D3q19, ViscousLattice3};
use crate::material::free_surface::sgs::SgsUnits;
use crate::material::free_surface::{
    MaterialFreeSurfaceSgsConfig3, MaterialFreeSurfaceSgsError3, MaterialFreeSurfaceSgsReport3,
};
use crate::sgs3::{SgsConfig, SgsError, SgsModel, smagorinsky_viscosity};

impl From<FreeSurfaceSgsError3> for Error {
    fn from(error: FreeSurfaceSgsError3) -> Self {
        Self::Sgs(error)
    }
}
impl From<MaterialFreeSurfaceSgsError3> for Error {
    fn from(error: MaterialFreeSurfaceSgsError3) -> Self {
        match error {
            MaterialFreeSurfaceSgsError3::Material(error) => error.into(),
            MaterialFreeSurfaceSgsError3::Sgs(error) => error.into(),
        }
    }
}

/// Source-resolved isothermal VOF with one immutable exploratory SGS closure.
///
/// Binding retains the original molecular law, physical scaling, pressure datum
/// and bulk selections. It prevents molecular samples fitted under total shear
/// relaxation from being consumed by a subsequent unbound molecular step.
/// No thermal deposition, subgrid interface model or turbulence qualification
/// is supplied. Configure the original physical bulk source before binding.
#[derive(Clone)]
pub struct MaterialFreeSurfaceRheologySgs3<'s, 'm> {
    flow: MaterialFreeSurfaceRheology3<'s, 'm>,
    config: MaterialFreeSurfaceSgsConfig3,
    native_config: SgsConfig,
}

impl<'s, 'm> MaterialFreeSurfaceRheology3<'s, 'm> {
    /// Bind SGS and admit the complete current source/total-relaxation field.
    /// Smagorinsky and disabled SGS need no gradient workspace; active WALE/Vreman require one with
    /// capacity for the current and intended future active tiles. The shared
    /// relaxation workspace is scratch, consumable only after success.
    /// The zero coefficient retains every original molecular source sample.
    pub fn with_sgs_si(
        self,
        config: MaterialFreeSurfaceSgsConfig3,
        gate: &CancelGate,
        kinematics: Option<&mut FreeSurfaceKinematicsWorkspace3>,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<MaterialFreeSurfaceRheologySgs3<'s, 'm>, Error> {
        poll(gate)?;
        let native_config = config.lattice(self.flow.scaling.cell_size_m)?;
        let mut result = MaterialFreeSurfaceRheologySgs3 {
            flow: self,
            config,
            native_config,
        };
        let (samples, _) =
            result.prepare(&result.flow.flow.lattice, gate, kinematics, workspace)?;
        if native_config.model.coefficient() != 0.0 {
            if let Some(field) = result.flow.flow.scaling.source_bulk_relaxations(
                samples
                    .iter()
                    .map(|(&location, cell)| (location, cell.material.material())),
                gate,
            )? {
                result
                    .flow
                    .flow
                    .lattice
                    .rebind_bulk_relaxations(field, gate)?;
            }
            result.flow.samples = samples;
        }
        poll(gate)?;
        Ok(result)
    }
}

impl MaterialFreeSurfaceRheologySgs3<'_, '_> {
    /// Original native owner, including the accepted force frame and bulk field.
    pub const fn lattice(&self) -> &FreeSurface3 {
        self.flow.lattice()
    }
    pub const fn config(&self) -> MaterialFreeSurfaceSgsConfig3 {
        self.config
    }
    pub fn material_at(&self, coordinates: [usize; 3]) -> Option<&ResolvedMaterialStatePoint> {
        self.flow.material_at(coordinates)
    }
    /// Accepted source shear rate and molecular tau. The SGS workspace exposes
    /// total collision tau after the successful endpoint preparation.
    pub fn rheology_at(&self, coordinates: [usize; 3]) -> Option<&IsothermalRheologySample> {
        self.flow.rheology_at(coordinates)
    }
    pub fn ledger_mass_kg(&self) -> Result<f64, Error> {
        self.flow.ledger_mass_kg()
    }
    pub fn sample_at_si(
        &self,
        x: usize,
        y: usize,
        z: usize,
    ) -> Result<Option<MaterialFreeSurfaceSample3>, Error> {
        self.flow.sample_at_si(x, y, z)
    }

    /// Advance one original VOF tick and re-admit the actual endpoint source.
    /// The returned powers describe the input state and use the constitutive
    /// stress-derived shear for all closures. WALE/Vreman's full gradient sets
    /// eddy viscosity, while its strain norm never replaces the molecular query.
    /// Newly wet cells resolve their own sources before bulk-field publication.
    pub fn step_with_gate(
        &mut self,
        gate: &CancelGate,
        kinematics: Option<&mut FreeSurfaceKinematicsWorkspace3>,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<MaterialFreeSurfaceSgsReport3, Error> {
        self.step_inner(gate, kinematics, workspace, None)
    }

    /// The same source/SGS step with positive prescribed absolute gas pressure
    /// [Pa] at exposed liquid centers [m]. The source uses solved liquid
    /// pressure and its original fixed datum, never the neighboring gas value.
    pub fn step_with_gas_pressure(
        &mut self,
        gate: &CancelGate,
        maximum_density_excursion: f64,
        pressure_pa: impl Fn([f64; 3]) -> f64,
        kinematics: Option<&mut FreeSurfaceKinematicsWorkspace3>,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<MaterialFreeSurfaceSgsReport3, Error> {
        poll(gate)?;
        let scale = self.flow.flow.scaling;
        let reference = scale.reference_pressure_pa.ok_or(Error::Invalid(
            "rheology SGS gas loading requires an explicit fixed pressure reference",
        ))?;
        let field = super::super::pressure::lattice_pressure_field(
            reference,
            scale.cell_size_m,
            scale.pressure_unit_pa,
            pressure_pa,
        );
        self.step_inner(
            gate,
            kinematics,
            workspace,
            Some((maximum_density_excursion, GasPressureField3::Cell(&field))),
        )
    }

    fn step_inner(
        &mut self,
        gate: &CancelGate,
        mut kinematics: Option<&mut FreeSurfaceKinematicsWorkspace3>,
        workspace: &mut FreeSurfaceSgsWorkspace3,
        gas_pressure: Option<(f64, GasPressureField3<'_>)>,
    ) -> Result<MaterialFreeSurfaceSgsReport3, Error> {
        poll(gate)?;
        let units = SgsUnits::new(&self.flow.flow.scaling)?;
        let (samples, report) =
            self.prepare(self.lattice(), gate, kinematics.as_deref_mut(), workspace)?;
        let report = units.report(report)?;
        if self.native_config.model.coefficient() == 0.0 {
            self.flow.step_inner(gate, gas_pressure)?;
            workspace.tau.retain(|&(key, lane), _| {
                self.flow
                    .samples
                    .contains_key(&MassLocation3::Cell(key, lane))
            });
            for (&location, cell) in &self.flow.samples {
                if let MassLocation3::Cell(key, lane) = location {
                    workspace.tau.insert((key, lane), cell.sample.collision_tau);
                }
            }
            return Ok(report);
        }
        let mut lattice = self.lattice().clone();
        self.rebind_bulk(&mut lattice, &samples, gate)?;
        drop(samples);
        lattice.step_with_material_relaxation_and_pressure(gate, &workspace.tau, gas_pressure)?;
        let (samples, _) = self.prepare(&lattice, gate, kinematics, workspace)?;
        self.rebind_bulk(&mut lattice, &samples, gate)?;
        poll(gate)?;
        self.flow.flow.lattice = lattice;
        self.flow.samples = samples;
        Ok(report)
    }

    fn rebind_bulk(
        &self,
        lattice: &mut FreeSurface3,
        samples: &BTreeMap<MassLocation3, CellMaterial>,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        if let Some(field) = self.flow.flow.scaling.source_bulk_relaxations(
            samples
                .iter()
                .map(|(&location, cell)| (location, cell.material.material())),
            gate,
        )? {
            lattice.rebind_bulk_relaxations(field, gate)?;
        }
        Ok(())
    }

    fn prepare(
        &self,
        lattice: &FreeSurface3,
        gate: &CancelGate,
        kinematics: Option<&mut FreeSurfaceKinematicsWorkspace3>,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<(BTreeMap<MassLocation3, CellMaterial>, FreeSurfaceSgsReport3), Error> {
        poll(gate)?;
        let scale = self.flow.flow.scaling;
        let config = self.native_config;
        molecular_relaxation(lattice.collision_model())?;
        let disabled = config.model.coefficient() == 0.0;
        let field =
            if disabled || matches!(config.model, SgsModel::Smagorinsky { .. }) {
                None
            } else {
                Some(lattice.kinematics(gate, kinematics.ok_or(Error::Invalid(
                "WALE/Vreman rheology requires the native liquid-gradient workspace",
            ))?).map_err(FreeSurfaceSgsError3::from)?)
            };
        let mass = lattice.mass_snapshot();
        workspace
            .tau
            .retain(|&(key, lane), _| mass.contains_key(&MassLocation3::Cell(key, lane)));
        let mut samples = BTreeMap::new();
        let mut report = FreeSurfaceSgsReport3::default();
        for (location, _) in mass {
            poll(gate)?;
            let MassLocation3::Cell(key, lane) = location else {
                continue;
            };
            let coordinates = FreeSurface3::transport_coordinates(location).expect("cell");
            let [x, y, z] = coordinates;
            let native = lattice
                .sample_at(x, y, z)
                .ok_or(Error::Invalid("missing wet SGS state"))?;
            let speed = native.velocity[0]
                .hypot(native.velocity[1])
                .hypot(native.velocity[2]);
            if !speed.is_finite() || speed > scale.flow.u_lattice {
                return Err(Error::Invalid(
                    "rheology SGS exceeds the admitted scenario speed",
                ));
            }
            let pressure = if self.flow.law.pressure_axis().is_some() {
                Some(scale.absolute_pressure_si(native.density)?)
            } else {
                None
            };
            let raw = lattice.transport_stress_numerators(location)?;
            let fixed_eddy = if let Some(field) = &field {
                if native.fill_fraction == 0.0 {
                    0.0
                } else {
                    gradient_sgs_at(field, coordinates, config)?.1
                }
            } else {
                0.0
            };
            let eddy_at = |rate: f64| -> Result<f64, Error> {
                Ok(match config.model {
                    SgsModel::Smagorinsky { cs } => {
                        smagorinsky_viscosity(cs, config.filter_width_lu, rate * scale.timestep_s)
                            .map_err(FreeSurfaceSgsError3::from)?
                    }
                    _ => fixed_eddy,
                })
            };
            let old = self.flow.samples.get(&location);
            let (cell, eddy, total_tau) = if disabled {
                let cell = old.ok_or(Error::Invalid(
                    "disabled SGS requires the accepted molecular sample",
                ))?;
                (Cow::Borrowed(cell), 0.0, cell.sample.collision_tau)
            } else {
                let initial = old.map_or(scale.flow.tau, |cell| cell.sample.collision_tau);
                let rate = old.map_or(raw[0] / initial / scale.timestep_s, |cell| {
                    cell.sample.shear_rate_per_s
                });
                let initial = initial + 3.0 * eddy_at(rate)?;
                let root = solve_local(
                    raw[0],
                    scale.timestep_s,
                    initial,
                    self.flow.solve,
                    || poll(gate),
                    |rate| {
                        let material = self.flow.law.resolve(rate, pressure)?;
                        scale.bulk_relaxation_rate_at(material.material())?;
                        let (plan, _) = plan_momentum_scaling(
                            material.kinematic_viscosity_m2_s(),
                            scale.cell_size_m,
                            scale.timestep_s,
                            0.0,
                        )?;
                        let eddy = eddy_at(rate)?;
                        let total = plan.tau + 3.0 * eddy;
                        if !total.is_finite() || total <= 0.5 {
                            return Err(Error::Invalid(
                                "unrepresentable molecular plus SGS relaxation",
                            ));
                        }
                        Ok(((material, plan.tau, eddy), total))
                    },
                )
                .map_err(|error| match error {
                    LocalRheologyError::Evaluation(error) => error,
                    LocalRheologyError::Invalid(what) => Error::Invalid(what),
                    LocalRheologyError::Convergence {
                        iterations,
                        relaxation_defect,
                    } => Error::Convergence {
                        coordinates,
                        iterations,
                        relaxation_defect,
                    },
                })?;
                let (material, molecular_tau, eddy) = root.value;
                (
                    Cow::Owned(CellMaterial {
                        material,
                        sample: IsothermalRheologySample {
                            shear_rate_per_s: root.shear_rate_per_s,
                            collision_tau: molecular_tau,
                            iterations: root.iterations,
                            relaxation_defect: root.relaxation_defect,
                        },
                    }),
                    eddy,
                    molecular_tau + 3.0 * eddy,
                )
            };
            let molecular_nu = (cell.sample.collision_tau - 0.5) / 3.0;
            let ratio = eddy / molecular_nu;
            if !ratio.is_finite() || ratio > config.maximum_eddy_to_molecular {
                return Err(
                    FreeSurfaceSgsError3::from(SgsError::ViscosityRatioExceeded {
                        cell: report.wet_cells,
                        ratio,
                        maximum: config.maximum_eddy_to_molecular,
                    })
                    .into(),
                );
            }
            D3q19::with_shear_relaxation(lattice.collision_model(), total_tau).map_err(
                |source| FreeSurfaceError3::Collision {
                    tile_key: key,
                    lane,
                    source,
                },
            )?;
            let bulk_rate = scale.bulk_relaxation_rate_at(cell.material.material())?;
            let bulk_nu = bulk_rate.map_or((2.0 / 3.0) * molecular_nu, |rate| {
                (2.0 / 9.0) * (1.0 / rate - 0.5)
            });
            let divergence = raw[1] * bulk_rate.unwrap_or(1.0 / total_tau) * (1.5_f64).sqrt();
            let rate = cell.sample.shear_rate_per_s * scale.timestep_s;
            let weight = native.fill_fraction * native.density;
            report.wet_cells += 1;
            report.liquid_volume_lu += native.fill_fraction;
            report.resolved_dissipation_lu += weight * molecular_nu * rate * rate;
            report.resolved_bulk_dissipation_lu += weight * bulk_nu * divergence * divergence;
            report.sgs_dissipation_lu += weight * eddy * rate * rate;
            report.maximum_eddy_viscosity_lu = report.maximum_eddy_viscosity_lu.max(eddy);
            report.maximum_eddy_to_molecular = report.maximum_eddy_to_molecular.max(ratio);
            workspace.tau.insert((key, lane), total_tau);
            if !disabled {
                samples.insert(location, cell.into_owned());
            }
        }
        if !report.resolved_dissipation_lu.is_finite()
            || !report.resolved_bulk_dissipation_lu.is_finite()
            || !report.sgs_dissipation_lu.is_finite()
        {
            return Err(Error::Invalid("unrepresentable rheology SGS dissipation"));
        }
        poll(gate)?;
        Ok((samples, report))
    }
}
