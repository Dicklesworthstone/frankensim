//! Local linear acoustic characteristics on an open NSF axis.
//!
//! For the outward normal, J+ = du_n + dp/(rho_ref c_ref) leaves the
//! domain and J- = du_n - dp/(rho_ref c_ref) enters it. Only J- is
//! prescribed. The entropy and tangential characteristics use their upwind
//! donor. This is the local linear, plane-normal rung of characteristic
//! boundaries (Poinsot & Lele, JCP 101 (1992), 104-129,
//! doi:10.1016/0021-9991(92)90046-2), not the full nonlinear NSCBC scheme.
//!
//! One reconstructed physical face state supplies the incoming lattice
//! populations, thermal-pressure force and total-energy flux. The incoming
//! populations extrapolate the adjacent post-collision non-equilibrium part;
//! no cell is overwritten with a reservoir and no sponge absorbs its energy.
//! A pressure-driven face instead prescribes dp and evaluates
//! J- = J+ - 2 dp/(rho_ref c_ref) at each original quadrature point and predictor.
//! The outgoing wave then determines normal velocity; the same reconstructed
//! face still owns populations, material response and conservative fluxes.
//! The normal boundary stencil is first order. Oblique transparency, a
//! <-40 dB reflection bound and resolved instrument accuracy are unqualified.

use super::{NsfError, NsfGrid3, NsfSgsTransport, Q, THETA0, check_state, equilibrium};
use crate::lattice3::Lattice3;

mod patch;
pub use patch::NsfCharacteristicPatch;
pub(super) use patch::CharacteristicPatchBinding;

#[cfg(test)]
mod thermal_tests;

/// Source-bound reference and one incoming plane-normal acoustic characteristic.
/// Both faces use their own outward normal; low is negative, high positive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfCharacteristicBoundary {
    pub reference_pressure_pa: f64,
    pub reference_temperature_k: f64,
    /// Cartesian velocity of the uniform exterior reference [m/s].
    pub reference_velocity_m_s: [f64; 3],
    /// `delta u_n - delta p / (rho_ref*c_ref)` [m/s]. Zero is the
    /// unperturbed incoming acoustic characteristic, not zero throughflow.
    /// This value is used by the incoming-characteristic drive. A prescribed
    /// pressure drive computes its local incoming value from the outgoing wave.
    pub incoming_characteristic_m_s: f64,
}

/// Physical control at one original characteristic face.
///
/// Pressure is relative to that face's retained material/reference pressure.
/// It is a uniform prescribed pressure, not a finite reservoir's dynamics,
/// stagnation pressure or imposed velocity. Explicit inflow temperature supplies
/// the incoming entropy characteristic through the same retained gas source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NsfBoundaryDrive {
    /// Hold J- = delta u_n - delta p/(rho_ref*c_ref) [m/s].
    IncomingCharacteristic(f64),
    /// Hold the reconstructed face pressure minus its reference pressure [Pa].
    /// Each local outgoing characteristic determines the resulting flow.
    PressureOffset(f64),
    /// Hold pressure offset [Pa] and incoming gas temperature [K]. Temperature
    /// is imposed only where the reconstructed normal flow is inward; outgoing
    /// and stationary faces retain the adjacent fluid's entropy characteristic.
    /// This is an ideal thermal source of the original fixed composition, not
    /// stagnation temperature, a finite player inventory or species injection.
    PressureTemperature {
        pressure_offset_pa: f64,
        inflow_temperature_k: f64,
    },
}

impl NsfBoundaryDrive {
    /// Scalar finiteness only; original face/material admission still applies.
    pub const fn is_finite(self) -> bool {
        match self {
            Self::IncomingCharacteristic(value) | Self::PressureOffset(value) => value.is_finite(),
            Self::PressureTemperature { pressure_offset_pa, inflow_temperature_k } =>
                pressure_offset_pa.is_finite() && inflow_temperature_k.is_finite(),
        }
    }
}

/// Two open faces on one axis of the original NSF grid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfCharacteristicAxis {
    pub axis: usize,
    /// Low face followed by high face.
    pub faces: [NsfCharacteristicBoundary; 2],
    /// Bound on |delta p|/p_ref, |delta T|/T_ref and |delta u|/c_ref
    /// at both the adjacent cells and reconstructed faces, in (0, 0.1].
    /// This admission bound limits the linearization domain; it is not an
    /// error estimate or a measured reflection guarantee.
    pub maximum_relative_perturbation: f64,
}

/// One accepted boundary quadrature cell, suitable for an explicit mortar.
/// Pressure and velocity are the actual common boundary state of the last
/// accepted predictor, including the original solver's finite coupling error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfBoundaryTrace {
    pub cell: [usize; 3],
    pub area_m2: f64,
    pub pressure_pa: f64,
    pub temperature_k: f64,
    pub density_kg_m3: f64,
    pub outward_normal_velocity_m_s: f64,
}

/// Exchange during one accepted lattice step. Positive means OUT of the grid.
/// The actual lattice mass flux is intentionally separate from the continuum
/// characteristic volume flow; their discretization defect is not hidden.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NsfBoundaryExchange {
    pub mass_out_kg: f64,
    /// Streamed momentum plus the trapezoidal thermal-pressure impulse.
    /// This includes pressure traction, not just advected momentum.
    pub momentum_out_kg_m_s: [f64; 3],
    /// Original FV enthalpy, conduction and viscous-work flux in the fixed
    /// sensible-energy chart of this grid's retained gas owner.
    pub total_energy_out_j: f64,
    /// Time-and-area averaged pressure of this face.
    pub mean_pressure_pa: f64,
    /// Time-averaged integral of the reconstructed normal velocity [m^3/s].
    pub normal_volume_flow_m3_s: f64,
    /// Integral of p*u_n over area and this step [J]. This is not generally
    /// mean_pressure_pa * normal_volume_flow_m3_s * dt for a nonuniform face.
    pub pressure_work_out_j: f64,
    /// Mean pressure relative to this face's declared reference [Pa].
    /// Subtraction occurs at the original quadrature, before time averaging.
    pub mean_acoustic_pressure_pa: f64,
    /// Integral of (p-p_ref)*u_n over the original face/time quadrature [J].
    /// Kept independently of the much larger absolute pressure work; it is
    /// not reconstructed by subtracting two nearly equal interval totals.
    pub acoustic_pressure_work_out_j: f64,
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct BoundaryFlux {
    pub total_out: [f64; 2],
    pub internal_out: [f64; 2],
    pub prescribed_out: f64,
    pub patch_total_out: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct OpenFace {
    pub rho: f64,
    pub theta: f64,
    pub u: [f64; 3],
    internal: f64,
    pressure_pa: f64,
}

#[derive(Clone, Debug)]
pub(super) struct OpenFaces {
    axis: usize,
    values: [Vec<OpenFace>; 2],
    patch: Option<CharacteristicPatchBinding>,
}

pub(super) struct IncomingPopulations {
    axis: usize,
    values: [Vec<[f64; Q]>; 2],
}

#[derive(Clone, Copy, Default)]
pub(super) struct LinkExchange {
    mass: f64,
    momentum: [f64; 3],
}

fn plane_index(x: [usize; 3], dims: [usize; 3], axis: usize) -> usize {
    match axis {
        0 => x[1] + dims[1] * x[2],
        1 => x[0] + dims[0] * x[2],
        _ => x[0] + dims[0] * x[1],
    }
}

impl IncomingPopulations {
    pub(super) fn at(&self, x: [usize; 3], dims: [usize; 3], side: usize) -> &[f64; Q] {
        &self.values[side][plane_index(x, dims, self.axis)]
    }
}

fn norm(u: [f64; 3]) -> f64 {
    (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt()
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Open both faces of one non-periodic axis before the first step.
    /// Transverse axes may be periodic or stationary halfway bounce-back walls.
    /// Two complete layers next to each open face must be fluid; interior
    /// stationary obstacles remain allowed. A diagonal link crossing both an
    /// open face and a transverse wall belongs to the wall, on every axis.
    /// Streaming and the open mass/momentum ledger use that same ownership.
    /// No extra cells, damping or second flow owner are introduced.
    ///
    /// Reference and face thermodynamics/transport query the SAME gas source
    /// as the interior. Frozen composition and the original source support,
    /// temperature, viscosity and Mach refusals all still apply.
    pub fn with_characteristic_axis(
        mut self,
        open: NsfCharacteristicAxis,
    ) -> Result<Self, NsfError> {
        if self.reservoir_faces.iter().flatten().any(Option::is_some) {
            return Err(NsfError::Input(
                "characteristic and prescribed reservoir faces are exclusive",
            ));
        }
        if self.steps != 0 || self.open.is_some() || self.sources.is_some() {
            return Err(NsfError::Input(
                "characteristic topology is cold-only and excludes manufactured body sources",
            ));
        }
        let axis = open.axis;
        if axis >= 3 || self.spec.periodic[axis] || self.spec.dims[axis] < 4 || self.moving_walls {
            return Err(NsfError::Input(
                "one non-periodic open axis needs >=4 cells and stationary transverse boundaries",
            ));
        }
        for i in 0..self.n {
            let x = self.at(i)[axis];
            if (x < 2 || x >= self.spec.dims[axis] - 2) && self.solid[i] {
                return Err(NsfError::Input(
                    "open faces require a two-cell all-fluid collar",
                ));
            }
        }
        self.rebind_characteristics(open, [None; 2], [None; 2], None)?;
        Ok(self)
    }

    /// Change only the incoming acoustic data at an accepted time boundary.
    /// Physical density, velocity and total energy are preserved. The Guo
    /// half-force population chart and current face-flux history are rebased
    /// together, so the next trapezoid uses the new data on both time levels.
    /// Failure publishes neither data nor partial boundary state.
    pub fn set_incoming_characteristics(&mut self, incoming_m_s: [f64; 2]) -> Result<(), NsfError> {
        self.set_boundary_drives(incoming_m_s.map(NsfBoundaryDrive::IncomingCharacteristic))
    }

    /// Change the physical drive at both original characteristic faces.
    ///
    /// A pressure drive fixes each reconstructed quadrature pressure throughout
    /// the original predictor/corrector. It does not freeze the incoming wave:
    /// the live outgoing wave determines its value and the normal flow. Inflow
    /// tangential data retain the original material reference. Incoming entropy
    /// retains that reference unless an explicit inflow temperature is supplied;
    /// outflow retains the adjacent cell. Pressure, temperature, Mach and source-domain
    /// refusals apply before any physical state or control is published.
    /// An attached aperture retains its own drive; these controls apply to the
    /// complementary cells of its face and every cell of the opposite face.
    pub fn set_boundary_drives(&mut self, drives: [NsfBoundaryDrive; 2]) -> Result<(), NsfError> {
        self.set_boundary_drives_with_binding(drives, self.characteristic_patch)
    }

    fn set_boundary_drives_with_binding(
        &mut self, drives: [NsfBoundaryDrive; 2], patch: Option<CharacteristicPatchBinding>,
    ) -> Result<(), NsfError> {
        let mut open = self
            .open
            .ok_or(NsfError::Input("no characteristic axis is attached"))?;
        let mut pressure = [None; 2];
        let mut temperature = [None; 2];
        for (side, drive) in drives.into_iter().enumerate() {
            if !drive.is_finite() {
                return Err(NsfError::Input("characteristic boundary drive must be finite"));
            }
            match drive {
                NsfBoundaryDrive::IncomingCharacteristic(value) => {
                    open.faces[side].incoming_characteristic_m_s = value;
                }
                NsfBoundaryDrive::PressureOffset(value) => {
                    // No single incoming characteristic represents a
                    // pressure-controlled face with a nonuniform outgoing wave.
                    open.faces[side].incoming_characteristic_m_s = 0.0;
                    pressure[side] = Some(value);
                }
                NsfBoundaryDrive::PressureTemperature { pressure_offset_pa, inflow_temperature_k } => {
                    open.faces[side].incoming_characteristic_m_s = 0.0;
                    pressure[side] = Some(pressure_offset_pa);
                    temperature[side] = Some(inflow_temperature_k);
                }
            }
        }
        self.rebind_characteristics(open, pressure, temperature, patch)
    }

    /// Actual accepted controls, low then high, with their physical units.
    /// An attached aperture has its separate `characteristic_patch_drive`.
    #[must_use]
    pub fn boundary_drives(&self) -> Option<[NsfBoundaryDrive; 2]> {
        let open = self.open.as_ref()?;
        Some(std::array::from_fn(|side| {
            self.characteristic_pressure_offset_pa[side].map_or(
                NsfBoundaryDrive::IncomingCharacteristic(
                    open.faces[side].incoming_characteristic_m_s,
                ),
                |pressure_offset_pa| self.characteristic_inflow_temperature_k[side].map_or(
                    NsfBoundaryDrive::PressureOffset(pressure_offset_pa),
                    |inflow_temperature_k| NsfBoundaryDrive::PressureTemperature {
                        pressure_offset_pa, inflow_temperature_k,
                    },
                ),
            )
        }))
    }

    #[must_use]
    pub const fn characteristic_axis(&self) -> Option<&NsfCharacteristicAxis> {
        self.open.as_ref()
    }

    /// The two face exchanges of the last accepted step, low then high.
    #[must_use]
    pub const fn last_boundary_exchange(&self) -> Option<&[NsfBoundaryExchange; 2]> {
        self.last_boundary_exchange.as_ref()
    }

    /// Quadrature traces in increasing original cell-index order.
    #[must_use]
    pub fn open_boundary_traces(&self, high: bool) -> Option<Vec<NsfBoundaryTrace>> {
        let faces = self.open_faces.as_ref()?;
        let side = usize::from(high);
        let sign = if high { 1.0 } else { -1.0 };
        let edge = if high {
            self.spec.dims[faces.axis] - 1
        } else {
            0
        };
        let mut traces = Vec::with_capacity(faces.values[side].len());
        for i in 0..self.n {
            let x = self.at(i);
            if x[faces.axis] != edge {
                continue;
            }
            let face = faces.values[side][plane_index(x, self.spec.dims, faces.axis)];
            traces.push(NsfBoundaryTrace {
                cell: x,
                area_m2: self.law.dx_m * self.law.dx_m,
                pressure_pa: face.pressure_pa,
                temperature_k: self.law.kelvin_of_theta(face.theta),
                density_kg_m3: face.rho * self.law.rho_ref_kg_m3,
                outward_normal_velocity_m_s: sign * face.u[faces.axis] * self.law.dx_m
                    / self.law.dt_s,
            });
        }
        Some(traces)
    }

    /// Complete original material gas at a declared characteristic reference.
    /// No source law, composition or transport coefficient is substituted.
    pub fn characteristic_reference_gas(
        &self,
        high: bool,
    ) -> Result<fs_material::gas::GasState, NsfError> {
        let open = self
            .open
            .ok_or(NsfError::Input("no characteristic axis is attached"))?;
        let face = open.faces[usize::from(high)];
        self.spec
            .gas
            .at(face.reference_temperature_k, face.reference_pressure_pa)
    }

    /// Actual allocated fluid/solid cell count, for a caller's trial budget.
    pub const fn cell_count(&self) -> usize {
        self.n
    }

    /// Fluid mass [kg], independent of the face bookkeeping.
    #[must_use]
    pub fn total_mass_kg(&self) -> f64 {
        let (rho, _) = self.macro_fields(&self.f);
        rho.iter()
            .zip(&self.solid)
            .filter_map(|(&r, &solid)| (!solid).then_some(r))
            .sum::<f64>()
            * self.law.rho_ref_kg_m3
            * self.law.dx_m.powi(3)
    }

    /// Physical, half-force-corrected fluid momentum [kg m/s].
    #[must_use]
    pub fn total_momentum_kg_m_s(&self) -> [f64; 3] {
        let (_, mom) = self.macro_fields(&self.f);
        let mut total = [0.0; 3];
        for i in 0..self.n {
            if self.solid[i] {
                continue;
            }
            for a in 0..3 {
                total[a] += mom[i][a] + 0.5 * self.force[i][a];
            }
        }
        let unit = self.law.rho_ref_kg_m3 * self.law.dx_m.powi(4) / self.law.dt_s;
        total.map(|p| p * unit)
    }

    fn rebind_characteristics(
        &mut self, open: NsfCharacteristicAxis, pressure_offset_pa: [Option<f64>; 2],
        inflow_temperature_k: [Option<f64>; 2],
        patch: Option<CharacteristicPatchBinding>,
    ) -> Result<(), NsfError> {
        let bound = open.maximum_relative_perturbation;
        if !bound.is_finite() || bound <= 0.0 || bound > 0.1 {
            return Err(NsfError::Input(
                "linear characteristic perturbation bound must lie in (0, 0.1]",
            ));
        }
        let (rho, mom) = self.macro_fields(&self.f);
        let u = self.velocity_field(&rho, &mom, &self.force);
        let faces = self.resolve_open_faces(Some(&open), pressure_offset_pa,
            inflow_temperature_k, patch, &rho, &u, &self.theta)?;
        let force = self.force_field(&rho, &self.theta, faces.as_ref());
        let mut f = self.f.clone();
        for i in 0..self.n {
            if self.solid[i] {
                continue;
            }
            super::rebase_force::<L>(
                &mut f[i], u[i],
                std::array::from_fn(|a| self.force[i][a] - force[i][a]),
            );
        }
        let transport = self.sgs_field(&rho, &u, &self.theta, &self.material, faces.as_ref())?;
        let sig = self.stress_field(&f, &force, &rho, &self.theta, &self.material, &transport);
        let (div, wall, flux) =
            self.flux_divergence(&self.e, &rho, &u, &self.theta, &sig, &transport, faces.as_ref())?;
        self.open = Some(open);
        self.characteristic_pressure_offset_pa = pressure_offset_pa;
        self.characteristic_inflow_temperature_k = inflow_temperature_k;
        self.characteristic_patch = patch;
        self.open_faces = faces;
        self.f = f;
        self.force = force;
        self.sgs_transport = transport;
        self.div = div;
        self.wall_out = wall;
        self.open_out = flux;
        self.last_boundary_exchange = None;
        self.last_patch_exchange = None;
        Ok(())
    }

    pub(super) fn resolve_open_faces(
        &self,
        open: Option<&NsfCharacteristicAxis>,
        pressure_offset_pa: [Option<f64>; 2],
        inflow_temperature_k: [Option<f64>; 2],
        patch: Option<CharacteristicPatchBinding>,
        rho: &[f64],
        u: &[[f64; 3]],
        theta: &[f64],
    ) -> Result<Option<OpenFaces>, NsfError> {
        let Some(open) = open else {
            return Ok(None);
        };
        let count = self.n / self.spec.dims[open.axis];
        let mut values = [Vec::with_capacity(count), Vec::with_capacity(count)];
        let velocity_unit = self.law.dx_m / self.law.dt_s;
        for (side, output) in values.iter_mut().enumerate() {
            let boundary = open.faces[side];
            let sign = if side == 0 { -1.0 } else { 1.0 };
            let edge = if side == 0 {
                0
            } else {
                self.spec.dims[open.axis] - 1
            };
            let reference = self.spec.gas.at(
                boundary.reference_temperature_k,
                boundary.reference_pressure_pa,
            )?;
            if !boundary.incoming_characteristic_m_s.is_finite()
                || !boundary
                    .reference_velocity_m_s
                    .iter()
                    .all(|v| v.is_finite())
                || boundary.reference_temperature_k < self.spec.temperature_window_k[0]
                || boundary.reference_temperature_k > self.spec.temperature_window_k[1]
                || norm(boundary.reference_velocity_m_s) / reference.sound_speed
                    > self.spec.maximum_mach
            {
                return Err(NsfError::Input(
                    "characteristic reference or incoming data outside the NSF domain",
                ));
            }
            let impedance = reference.density * reference.sound_speed;
            let alpha = (reference.gamma - 1.0) / reference.gamma * reference.temperature
                / reference.pressure;
            let normal_reference = sign * boundary.reference_velocity_m_s[open.axis];
            // Admit a declared thermal source even while this face is outgoing:
            // a reversal must not activate previously hidden invalid material.
            // Legacy controls perform no additional material query.
            let admit_temperature = |offset: f64, temperature: f64| -> Result<(), NsfError> {
                if !temperature.is_finite()
                    || temperature < self.spec.temperature_window_k[0]
                    || temperature > self.spec.temperature_window_k[1]
                    || (temperature / reference.temperature - 1.0).abs()
                        > open.maximum_relative_perturbation
                {
                    return Err(NsfError::Input(
                        "characteristic inflow temperature leaves its declared source domain",
                    ));
                }
                self.spec.gas.at(temperature, reference.pressure + offset)?;
                Ok(())
            };
            if let Some(temperature) = inflow_temperature_k[side] {
                admit_temperature(pressure_offset_pa[side].ok_or(NsfError::Input(
                    "characteristic inflow temperature requires an explicit pressure control",
                ))?, temperature)?;
            }
            if let Some(patch) = patch.filter(|patch| usize::from(patch.high) == side) {
                if let NsfBoundaryDrive::PressureTemperature {
                    pressure_offset_pa, inflow_temperature_k,
                } = patch.drive {
                    admit_temperature(pressure_offset_pa, inflow_temperature_k)?;
                }
            }
            for i in 0..self.n {
                let x = self.at(i);
                if x[open.axis] != edge {
                    continue;
                }
                let pi = self.law.pressure_pa(rho[i], theta[i]);
                let ti = self.law.kelvin_of_theta(theta[i]);
                let ui = u[i].map(|v| v * velocity_unit);
                let outgoing =
                    sign * ui[open.axis] - normal_reference + (pi - reference.pressure) / impedance;
                let (incoming, pressure_offset, inflow_temperature) = match patch.filter(|patch|
                    usize::from(patch.high) == side && patch.shape.contains(x, open.axis))
                {
                    Some(patch) => match patch.drive {
                        NsfBoundaryDrive::IncomingCharacteristic(value) => (value, None, None),
                        NsfBoundaryDrive::PressureOffset(value) => (0.0, Some(value), None),
                        NsfBoundaryDrive::PressureTemperature { pressure_offset_pa, inflow_temperature_k } =>
                            (0.0, Some(pressure_offset_pa), Some(inflow_temperature_k)),
                    },
                    None => (boundary.incoming_characteristic_m_s, pressure_offset_pa[side],
                        inflow_temperature_k[side]),
                };
                let (normal, pressure) = match pressure_offset {
                    Some(offset) => {
                        // J- = J+ - 2 dp/Z. Evaluate the equivalent normal
                        // velocity directly, avoiding cancellation between
                        // two similar characteristic values. The pressure is
                        // exact to addition roundoff at every original face.
                        (normal_reference + outgoing - offset / impedance,
                            reference.pressure + offset)
                    }
                    None => (
                        normal_reference + 0.5 * (outgoing + incoming),
                        reference.pressure
                            + 0.5 * impedance * (outgoing - incoming),
                    ),
                };
                let entropy_temperature = if normal >= 0.0 {
                    ti - reference.temperature - alpha * (pi - reference.pressure)
                } else {
                    0.0
                };
                let temperature = match inflow_temperature.filter(|_| normal < 0.0) {
                    Some(temperature) => temperature,
                    None => reference.temperature
                        + alpha * (pressure - reference.pressure)
                        + entropy_temperature,
                };
                let mut velocity = if normal >= 0.0 {
                    ui
                } else {
                    boundary.reference_velocity_m_s
                };
                velocity[open.axis] = sign * normal;
                let du_i = std::array::from_fn(|a| ui[a] - boundary.reference_velocity_m_s[a]);
                let du_f =
                    std::array::from_fn(|a| velocity[a] - boundary.reference_velocity_m_s[a]);
                let perturbation = ((pi / reference.pressure) - 1.0)
                    .abs()
                    .max((ti / reference.temperature - 1.0).abs())
                    .max(norm(du_i) / reference.sound_speed)
                    .max((pressure / reference.pressure - 1.0).abs())
                    .max((temperature / reference.temperature - 1.0).abs())
                    .max(norm(du_f) / reference.sound_speed);
                if !perturbation.is_finite() || perturbation > open.maximum_relative_perturbation {
                    return Err(NsfError::Input(
                        "characteristic cell or face leaves its declared linearization domain",
                    ));
                }
                let face_rho =
                    pressure / (self.law.gas_constant * temperature * self.law.rho_ref_kg_m3);
                let face_theta = self.law.theta_of_kelvin(temperature);
                let face_u = velocity.map(|v| v / velocity_unit);
                let material = self.spec.gas.resolve(&self.law, face_rho, face_theta)?;
                let gamma = match &material {
                    Some(m) => m.state.gamma,
                    None => self.spec.gas.at(temperature, pressure)?.gamma,
                };
                check_state(&self.spec, gamma, i, face_rho, face_u, face_theta)?;
                if temperature < self.spec.temperature_window_k[0]
                    || temperature > self.spec.temperature_window_k[1]
                {
                    return Err(NsfError::Window {
                        cell: i,
                        temperature_k: temperature,
                        low_k: self.spec.temperature_window_k[0],
                        high_k: self.spec.temperature_window_k[1],
                    });
                }
                let internal = material
                    .as_ref()
                    .map_or(face_theta / (self.law.gamma - 1.0), |m| {
                        m.sensible_energy_lattice
                    });
                output.push(OpenFace {
                    rho: face_rho,
                    theta: face_theta,
                    u: face_u,
                    internal,
                    pressure_pa: pressure,
                });
            }
        }
        Ok(Some(OpenFaces {
            axis: open.axis,
            values,
            patch,
        }))
    }

    pub(super) fn characteristic_face<'a>(
        &self,
        faces: Option<&'a OpenFaces>,
        i: usize,
        axis: usize,
        direction: i64,
    ) -> Option<&'a OpenFace> {
        let faces = faces?;
        if axis != faces.axis {
            return None;
        }
        let side = usize::from(direction > 0);
        let x = self.at(i);
        let edge = if side == 0 {
            0
        } else {
            self.spec.dims[axis] - 1
        };
        if x[axis] != edge {
            return None;
        }
        Some(&faces.values[side][plane_index(x, self.spec.dims, axis)])
    }

    pub(super) fn open_energy_flux(
        &self,
        face: &OpenFace,
        i: usize,
        axis: usize,
        direction: i64,
        rho: &[f64],
        theta: &[f64],
        sig: &[[f64; 6]],
        transport: &[NsfSgsTransport],
    ) -> Result<(f64, f64), NsfError> {
        let mut k = self.spec.gas.face_conductivity(
            &self.law,
            0.5 * (theta[i] + face.theta),
            0.5 * (self.law.pressure_pa(rho[i], theta[i]) + face.pressure_pa),
        )?;
        // The modeled coefficient is extrapolated from this adjacent cell;
        // the original characteristic temperature owns the half-cell gradient.
        if let Some(local) = transport.get(i) {
            k += local.turbulent_conductivity_lu;
        }
        let heat = -2.0 * k * direction as f64 * (face.theta - theta[i]);
        let s = sig[i];
        let row = match axis {
            0 => [s[0], s[3], s[4]],
            1 => [s[3], s[1], s[5]],
            _ => [s[4], s[5], s[2]],
        };
        let viscous_work = row[0] * face.u[0] + row[1] * face.u[1] + row[2] * face.u[2];
        let kinetic = 0.5 * (face.u[0] * face.u[0] + face.u[1] * face.u[1] + face.u[2] * face.u[2]);
        let total =
            face.rho * (face.internal + kinetic + face.theta) * face.u[axis] - viscous_work + heat;
        let internal = face.rho * face.internal * face.u[axis] + heat;
        Ok((total, internal))
    }

    pub(super) fn open_incoming_populations(
        &self,
        post: &[[f64; Q]],
        rho: &[f64],
        u: &[[f64; 3]],
    ) -> Option<IncomingPopulations> {
        let faces = self.open_faces.as_ref()?;
        let mut values = [
            Vec::with_capacity(faces.values[0].len()),
            Vec::with_capacity(faces.values[1].len()),
        ];
        for (side, output) in values.iter_mut().enumerate() {
            let edge = if side == 0 {
                0
            } else {
                self.spec.dims[faces.axis] - 1
            };
            for i in 0..self.n {
                let x = self.at(i);
                if x[faces.axis] != edge {
                    continue;
                }
                let face = faces.values[side][plane_index(x, self.spec.dims, faces.axis)];
                let target = equilibrium::<L>(face.rho, face.u, THETA0);
                let interior = equilibrium::<L>(rho[i], u[i], THETA0);
                output.push(std::array::from_fn(|q| {
                    target[q] + (post[i][q] - interior[q])
                }));
            }
        }
        Some(IncomingPopulations {
            axis: faces.axis,
            values,
        })
    }

    /// A link touching a stationary transverse wall is reflected there even
    /// when it also crosses the open plane. The inverse link has the same
    /// owner, so a corner can neither create an unmatched population nor leak
    /// through a wall. `direction` is +1 for outgoing and -1 for pull streaming.
    pub(super) fn characteristic_link_is_open(
        &self,
        x: [usize; 3],
        c: [i32; 3],
        direction: i32,
    ) -> bool {
        let Some(open) = self.open else {
            return false;
        };
        (0..3).all(|axis| {
            if axis == open.axis || self.spec.periodic[axis] {
                return true;
            }
            let next = x[axis] as i64 + i64::from(direction) * i64::from(c[axis]);
            next >= 0 && next < self.spec.dims[axis] as i64
        })
    }

    pub(super) fn open_link_exchange(
        &self,
        post: &[[f64; Q]],
        incoming: Option<&IncomingPopulations>,
    ) -> Option<([LinkExchange; 2], Option<LinkExchange>)> {
        let incoming = incoming?;
        let axis = incoming.axis;
        let mut result = [LinkExchange::default(); 2];
        let mut patch_exchange = self.characteristic_patch.map(|_| LinkExchange::default());
        for (side, output) in result.iter_mut().enumerate() {
            let sign = if side == 0 { -1 } else { 1 };
            let edge = if side == 0 {
                0
            } else {
                self.spec.dims[axis] - 1
            };
            for i in 0..self.n {
                let x = self.at(i);
                if x[axis] != edge {
                    continue;
                }
                let ghost = incoming.at(x, self.spec.dims, side);
                for q in 0..Q {
                    let c = L::E[q];
                    let normal = i32::from(c[axis]) * sign;
                    let exchanged = if normal > 0 && self.characteristic_link_is_open(x, c, 1) {
                        post[i][q]
                    } else if normal < 0 && self.characteristic_link_is_open(x, c, -1) {
                        -ghost[q]
                    } else {
                        continue;
                    };
                    output.mass += exchanged;
                    for a in 0..3 {
                        output.momentum[a] += f64::from(c[a]) * exchanged;
                    }
                    if self.characteristic_patch.is_some_and(|patch|
                        usize::from(patch.high) == side && patch.shape.contains(x, axis))
                    {
                        let patch = patch_exchange.as_mut().expect("selected patch");
                        patch.mass += exchanged;
                        for a in 0..3 {
                            patch.momentum[a] += f64::from(c[a]) * exchanged;
                        }
                    }
                }
            }
        }
        Some((result, patch_exchange))
    }

    pub(super) fn finish_boundary_exchange(
        &self,
        links: Option<([LinkExchange; 2], Option<LinkExchange>)>,
        before: Option<&OpenFaces>,
        after: Option<&OpenFaces>,
        flux_after: BoundaryFlux,
    ) -> Option<[NsfBoundaryExchange; 2]> {
        let (links, _) = links?;
        let (before, after) = (before?, after?);
        let axis = before.axis;
        let mass_unit = self.law.rho_ref_kg_m3 * self.law.dx_m.powi(3);
        let velocity_unit = self.law.dx_m / self.law.dt_s;
        let area = self.law.dx_m * self.law.dx_m;
        let mut result = [NsfBoundaryExchange::default(); 2];
        for side in 0..2 {
            let sign = if side == 0 { -1.0 } else { 1.0 };
            let out = &mut result[side];
            let reference_pressure = self.open.as_ref()?.faces[side].reference_pressure_pa;
            out.mass_out_kg = links[side].mass * mass_unit;
            out.momentum_out_kg_m_s = links[side].momentum.map(|p| p * mass_unit * velocity_unit);
            out.total_energy_out_j = 0.5
                * (self.open_out.total_out[side] + flux_after.total_out[side])
                * self.law.energy_unit_j();
            for (old, new) in before.values[side].iter().zip(&after.values[side]) {
                out.momentum_out_kg_m_s[axis] += 0.5
                    * sign
                    * (old.rho * (old.theta - THETA0) + new.rho * (new.theta - THETA0))
                    * mass_unit
                    * velocity_unit;
                out.mean_pressure_pa += 0.5 * (old.pressure_pa + new.pressure_pa);
                let old_acoustic = old.pressure_pa - reference_pressure;
                let new_acoustic = new.pressure_pa - reference_pressure;
                out.mean_acoustic_pressure_pa += 0.5 * (old_acoustic + new_acoustic);
                out.acoustic_pressure_work_out_j += 0.5
                    * sign
                    * (old_acoustic * old.u[axis] + new_acoustic * new.u[axis])
                    * velocity_unit
                    * area
                    * self.law.dt_s;
                out.normal_volume_flow_m3_s +=
                    0.5 * sign * (old.u[axis] + new.u[axis]) * velocity_unit * area;
                out.pressure_work_out_j += 0.5
                    * sign
                    * (old.pressure_pa * old.u[axis] + new.pressure_pa * new.u[axis])
                    * velocity_unit
                    * area
                    * self.law.dt_s;
            }
            out.mean_pressure_pa /= before.values[side].len() as f64;
            out.mean_acoustic_pressure_pa /= before.values[side].len() as f64;
        }
        Some(result)
    }
}
