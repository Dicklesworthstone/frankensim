//! Finite-source illumination of closed rigid or admittance scatterers.
//!
//! For exp(-i omega t), a monopole with outward volume velocity Q has
//! `p=-i*omega*rho*Q*exp(ikr)/(4*pi*r)` and
//! `v=Q*exp(ikr)/(4*pi*r)*(1/r-i*k)*r_hat`. The reactive 1/r² velocity
//! is retained: a nearby emitter is not a plane wave aimed at a centroid.
//! On a stationary wall, `v_total=-Y*p_total`, hence the scattered solve's
//! drive is `-v_incident-Y*p_incident`. Both original BM rows receive it.
//!
//! The source may be outside any number of admitted closed components.
//! Interior sources, surface-touching sources and singular receivers refuse.
//! This is homogeneous lossless scalar acoustics, not a finite-size radiator,
//! a fluid/solid transmission model or a calibrated anatomical HRTF.

use super::*;
use crate::near_field;
use std::f64::consts::PI;

#[path = "helmholtz_incident_traces.rs"]
mod traces;
pub use traces::IncidentTraces;

/// Maximum number of coherent point sources in one incident field.
pub const MAX_MONOPOLES: usize = near_field::MAX_RECEIVERS;

/// Maximum total number of coherent plane-wave and monopole source terms.
pub const MAX_INCIDENT_SOURCES: usize = MAX_MONOPOLES;

/// Ideal point emitter in the same stationary medium as the scatterer.
#[derive(Clone, Copy, Debug)]
pub struct Monopole {
    /// Physical source position [m].
    pub position_m: [f64; 3],
    /// Complex outward volume flow [m³/s], with exp(-i omega t).
    pub volume_velocity_m3_s: C64,
}

/// One analytic source term at the common incident frequency. Combining
/// plane waves and nearby emitters retains their complex interference.
#[derive(Clone, Copy, Debug)]
pub enum IncidentSource {
    PlaneWave(PlaneWave),
    Monopole(Monopole),
}

/// Coherent incident field. Sources are summed before scattering; their
/// phases and volume flows are not replaced by incoherent source power.
#[derive(Clone, Copy)]
pub enum IncidentField<'a> {
    /// Unbounded incident plane wave with pressure amplitude [Pa].
    PlaneWave(PlaneWave),
    /// Coherent finite emitters with individually specified volume flows.
    Monopoles(&'a [Monopole]),
    /// A coherent mixture of plane waves and finite-distance emitters.
    Sources(&'a [IncidentSource]),
}

/// Pressure and Euler particle velocity of the unscattered incident field.
#[derive(Clone, Copy, Debug)]
pub struct IncidentSample {
    /// Complex pressure [Pa], exp(-i omega t).
    pub pressure_pa: C64,
    /// Cartesian fluid particle velocity [m/s].
    pub velocity_m_s: [C64; 3],
}

pub(super) fn validate_frequency(k: f64, medium: Medium) -> Result<(), HelmholtzError> {
    let omega_rho = k * medium.sound_speed * medium.density;
    let rho_c = medium.sound_speed * medium.density;
    if !k.is_finite()
        || k <= 0.0
        || !medium.density.is_finite()
        || medium.density <= 0.0
        || !medium.sound_speed.is_finite()
        || medium.sound_speed <= 0.0
        || !omega_rho.is_finite()
        || omega_rho <= 0.0
        || !rho_c.is_finite()
        || rho_c <= 0.0
    {
        return Err(bad(
            "incident field requires finite positive k, density and sound speed",
        ));
    }
    Ok(())
}

impl IncidentField<'_> {
    fn validate(self) -> Result<(), HelmholtzError> {
        match self {
            Self::PlaneWave(plane) => {
                if plane.direction.iter().any(|v| !v.is_finite())
                    || (dot(plane.direction, plane.direction) - 1.0).abs() > 1e-10
                    || !finite(plane.amplitude_pa)
                {
                    return Err(bad(
                        "incident plane wave requires a finite amplitude and unit direction",
                    ));
                }
            }
            Self::Monopoles(sources) => {
                if sources.is_empty()
                    || sources.len() > MAX_MONOPOLES
                    || sources.iter().any(|s| {
                        s.position_m.iter().any(|v| !v.is_finite() || v.abs() > 1e9)
                            || !finite(s.volume_velocity_m3_s)
                    })
                {
                    return Err(bad(
                        "incident monopoles require a nonempty bounded finite source array",
                    ));
                }
            }
            Self::Sources(sources) => {
                if sources.is_empty() || sources.len() > MAX_INCIDENT_SOURCES {
                    return Err(bad("coherent source array must be nonempty and bounded"));
                }
                for source in sources {
                    match source {
                        IncidentSource::PlaneWave(plane) => Self::PlaneWave(*plane).validate()?,
                        IncidentSource::Monopole(source) => {
                            Self::Monopoles(std::slice::from_ref(source)).validate()?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Exact analytic incident pressure and velocity at a nonsingular point.
    /// This evaluates free space; exterior-to-scatterer admission belongs to
    /// `solve_scattering` and its batch form.
    pub fn pressure_velocity_at(
        self,
        k: f64,
        medium: Medium,
        position_m: [f64; 3],
    ) -> Result<IncidentSample, HelmholtzError> {
        validate_frequency(k, medium)?;
        self.validate()?;
        self.sample(k, medium, position_m)
    }

    fn sample(
        self,
        k: f64,
        medium: Medium,
        point: [f64; 3],
    ) -> Result<IncidentSample, HelmholtzError> {
        if point.iter().any(|v| !v.is_finite() || v.abs() > 1e9) {
            return Err(bad(
                "incident receiver must have finite resolved coordinates",
            ));
        }
        let mut out = IncidentSample {
            pressure_pa: C64::ZERO,
            velocity_m_s: [C64::ZERO; 3],
        };
        match self {
            Self::PlaneWave(plane) => {
                let phase = k * dot(plane.direction, point);
                if !phase.is_finite() {
                    return Err(bad("incident plane-wave phase overflow"));
                }
                out.pressure_pa = plane.amplitude_pa * expik(phase);
                out.velocity_m_s = plane.direction.map(|d| {
                    out.pressure_pa
                        .scale(d / (medium.density * medium.sound_speed))
                });
            }
            Self::Monopoles(sources) => {
                let omega_rho = k * medium.sound_speed * medium.density;
                for source in sources {
                    let delta = sub(point, source.position_m);
                    let radius = norm(delta);
                    let coordinate_scale = point
                        .into_iter()
                        .chain(source.position_m)
                        .fold(0.0_f64, |m, v| m.max(v.abs()));
                    if !radius.is_finite()
                        || radius <= 128.0 * f64::EPSILON * coordinate_scale
                        || radius <= 0.0
                        || !(k * radius).is_finite()
                    {
                        return Err(bad(
                            "monopole receiver coincides with its source or is unresolved",
                        ));
                    }
                    let flow_green = source.volume_velocity_m3_s
                        * expik(k * radius).scale(1.0 / (4.0 * PI * radius));
                    out.pressure_pa = out.pressure_pa + flow_green * C64::new(0.0, -omega_rho);
                    let radial = flow_green * C64::new(1.0 / radius, -k);
                    for (velocity, distance) in out.velocity_m_s.iter_mut().zip(delta) {
                        *velocity = *velocity + radial.scale(distance / radius);
                    }
                }
            }
            Self::Sources(sources) => {
                for source in sources {
                    let sample = match source {
                        IncidentSource::PlaneWave(plane) => {
                            Self::PlaneWave(*plane).sample(k, medium, point)?
                        }
                        IncidentSource::Monopole(source) => {
                            Self::Monopoles(std::slice::from_ref(source))
                                .sample(k, medium, point)?
                        }
                    };
                    out.pressure_pa = out.pressure_pa + sample.pressure_pa;
                    for (velocity, value) in out.velocity_m_s.iter_mut().zip(sample.velocity_m_s) {
                        *velocity = *velocity + value;
                    }
                }
            }
        }
        if !finite(out.pressure_pa) || out.velocity_m_s.iter().any(|v| !finite(*v)) {
            return Err(bad("incident pressure or velocity overflow"));
        }
        Ok(out)
    }

    fn boundary_drive(
        self,
        surface: &SpherePanels,
        k: f64,
        medium: Medium,
        law: Admittance<'_>,
        minimum_source_clearance_m: f64,
    ) -> Result<(Vec<C64>, Vec<C64>, Vec<C64>), HelmholtzError> {
        self.validate()?;
        if !minimum_source_clearance_m.is_finite() || minimum_source_clearance_m <= 0.0 {
            return Err(bad("incident source clearance must be positive and finite"));
        }
        if let Self::PlaneWave(plane) = self {
            // Preserve the established plane-wave arithmetic and regression bits.
            let (pressure, drive) = scattering_drive(surface, k, medium, law, plane)?;
            let velocity = plane_wave_normal_velocity(surface, medium, plane, &pressure);
            return Ok((pressure, velocity, drive));
        }
        let traces =
            IncidentTraces::from_field(surface, k, medium, self, minimum_source_clearance_m)?;
        let drive = traces.drive(law)?;
        let (pressure, velocity) = traces.into_boundary_traces();
        Ok((pressure, velocity, drive))
    }
}

/// Exact coherent free-space pressure and velocity of the declared sources.
/// Body-source clearance is admitted by [`IncidentTraces::from_sources`].
pub fn sample_incident_sources(
    k: f64,
    medium: Medium,
    sources: &[IncidentSource],
    point_m: [f64; 3],
) -> Result<IncidentSample, HelmholtzError> {
    IncidentField::Sources(sources).pressure_velocity_at(k, medium, point_m)
}

/// One solved scattering problem with its immutable incident-source and body
/// identity. The outgoing scattered trace is kept separate from total pressure.
pub struct ScatteredField<'s, 'i> {
    surface: &'s SpherePanels,
    incident: IncidentField<'i>,
    solution: ScatteringSolution,
}

impl ScatteredField<'_, '_> {
    /// Original outgoing boundary solution and total-wall absorption.
    pub fn solution(&self) -> &ScatteringSolution {
        &self.solution
    }

    /// Physical incident+diffracted pressure at finite exterior receivers.
    /// Quadrature estimates cover the scattered representation only, not its
    /// spatial discretization or the ideal point-source model.
    pub fn total_pressure_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_receiver_clearance_m: f64,
        options: near_field::Options,
    ) -> Result<near_field::Observation, HelmholtzError> {
        let trace = &self.solution.scattered;
        if points_m.is_empty() || points_m.len() > near_field::MAX_RECEIVERS {
            return Err(bad(
                "total receiver array exceeds the near-field cardinality limit",
            ));
        }
        // Refuse source collisions before the more expensive BEM quadrature.
        let incident: Vec<_> = points_m
            .iter()
            .map(|&x| self.incident.sample(trace.k, trace.medium, x))
            .collect::<Result<_, _>>()?;
        let geometry =
            near_field::Geometry::new(self.surface, points_m, minimum_receiver_clearance_m)?;
        let mut observation = geometry
            .prepare(trace.k, trace.medium, options)?
            .evaluate(trace)?;
        for (p, source) in observation.pressure.iter_mut().zip(incident) {
            *p = *p + source.pressure_pa;
            if !finite(*p) {
                return Err(bad("total scattering receiver pressure overflow"));
            }
        }
        Ok(observation)
    }

    /// Outgoing incident+scattered amplitude F [Pa*m], with p~F exp(ikr)/r.
    /// A plane wave has no outgoing 1/r total far field and is refused here.
    pub fn total_far_field(&self, direction: [f64; 3]) -> Result<C64, HelmholtzError> {
        if matches!(self.incident, IncidentField::PlaneWave(_))
            || matches!(self.incident, IncidentField::Sources(sources)
                if sources.iter().any(|s| matches!(s, IncidentSource::PlaneWave(_))))
        {
            return Err(bad(
                "a plane wave has no outgoing total far-field amplitude",
            ));
        }
        if direction.iter().any(|v| !v.is_finite()) || (norm(direction) - 1.0).abs() > 1e-10 {
            return Err(bad("total far field requires a unit direction"));
        }
        let trace = &self.solution.scattered;
        let mut value = far_field_integrated(self.surface, trace, trace.medium, &[direction])[0];
        let factor = C64::new(
            0.0,
            -trace.k * trace.medium.sound_speed * trace.medium.density / (4.0 * PI),
        );
        let mut add_source = |source: &Monopole| {
            value = value
                + factor
                    * source.volume_velocity_m3_s
                    * expik(-trace.k * dot(direction, source.position_m));
        };
        match self.incident {
            IncidentField::Monopoles(sources) => sources.iter().for_each(&mut add_source),
            IncidentField::Sources(sources) => {
                for source in sources {
                    if let IncidentSource::Monopole(source) = source {
                        add_source(source);
                    }
                }
            }
            IncidentField::PlaneWave(_) => unreachable!("plane waves refused above"),
        }
        if !finite(value) {
            return Err(bad("total far-field amplitude overflow"));
        }
        Ok(value)
    }
}

/// Solve exact finite-source illumination through the existing dense Robin
/// owner. The body may include disconnected scatterers; their interaction is
/// retained by one combined boundary solve, not by independent shadow factors.
pub fn solve_scattering<'s, 'i>(
    surface: &'s SpherePanels,
    k: f64,
    medium: Medium,
    law: Admittance<'_>,
    incident: IncidentField<'i>,
    minimum_source_clearance_m: f64,
    formulation: Formulation,
) -> Result<ScatteredField<'s, 'i>, HelmholtzError> {
    Ok(solve_scattering_batch(
        surface,
        k,
        medium,
        law,
        &[incident],
        minimum_source_clearance_m,
        formulation,
    )?
    .remove(0))
}

impl<'s> PreparedBoundary<'s> {
    /// Scatter a new analytic field using this existing factorization, while
    /// retaining its source identity for coherent finite receiver observations.
    pub fn scatter_field<'i>(
        &self,
        incident: IncidentField<'i>,
        minimum_source_clearance_m: f64,
    ) -> Result<ScatteredField<'s, 'i>, HelmholtzError> {
        let (pressure, velocity, drive) = incident.boundary_drive(
            self.surface,
            self.k,
            self.medium,
            self.law,
            minimum_source_clearance_m,
        )?;
        let solution = self.solve(&drive)?;
        Ok(ScatteredField {
            surface: self.surface,
            incident,
            solution: finish_scattering(self.surface, solution, pressure, velocity, self.law)?,
        })
    }
}

/// Finite-source or mixed illumination through the existing FMM/GMRES owner.
/// Source admission supports its original surface range; finite receivers
/// retain their independent dense observer cap. The original convergence,
/// cancellation and work receipts remain authoritative.
#[allow(clippy::too_many_arguments)]
pub fn solve_scattering_iterative<'s, 'i>(
    surface: &'s SpherePanels,
    k: f64,
    medium: Medium,
    law: Admittance<'_>,
    incident: IncidentField<'i>,
    minimum_source_clearance_m: f64,
    formulation: Formulation,
    options: &IterativeOptions,
    gate: &CancelGate,
) -> Result<
    (
        ScatteredField<'s, 'i>,
        fs_solver::krylov::SolveReport,
        ResourceReceipt,
    ),
    IterativeError,
> {
    if gate.is_requested() {
        return Err(IterativeError::Cancelled { matvecs: 0 });
    }
    law.validate(surface.areas().len())?;
    let traces =
        IncidentTraces::from_field(surface, k, medium, incident, minimum_source_clearance_m)?;
    let (solution, report, receipt) =
        solve_incident_bc_iterative(surface, k, medium, law, &traces, formulation, options, gate)?;
    if gate.is_requested() {
        return Err(IterativeError::Cancelled {
            matvecs: receipt.matvecs,
        });
    }
    Ok((
        ScatteredField {
            surface,
            incident,
            solution,
        },
        report,
        receipt,
    ))
}

/// Independent incident fields through one boundary assembly and LU. Source
/// admission precedes assembly; all RHS pass the original boundary work checks.
/// Dense assembly/LU inherit the original owner's cancellation limitation.
pub fn solve_scattering_batch<'s, 'i>(
    surface: &'s SpherePanels,
    k: f64,
    medium: Medium,
    law: Admittance<'_>,
    incidents: &[IncidentField<'i>],
    minimum_source_clearance_m: f64,
    formulation: Formulation,
) -> Result<Vec<ScatteredField<'s, 'i>>, HelmholtzError> {
    validate_frequency(k, medium)?;
    if incidents.is_empty() || incidents.len() > MAX_RADIATION_FIELDS_PER_BATCH {
        return Err(bad(
            "incident field batch must be nonempty and within the radiation batch cap",
        ));
    }
    admit_resolution(surface, k)?;
    validate_closed(surface)?;
    law.validate(surface.areas().len())?;
    let traces: Vec<_> = incidents
        .iter()
        .map(|&field| field.boundary_drive(surface, k, medium, law, minimum_source_clearance_m))
        .collect::<Result<_, _>>()?;
    let drives: Vec<_> = traces.iter().map(|(_, _, v)| v.as_slice()).collect();
    let solved = solve_radiation_bc_batch(surface, k, medium, &drives, law, formulation)?;
    solved
        .into_iter()
        .zip(traces)
        .zip(incidents)
        .map(|((solution, (p, v, _)), &incident)| {
            Ok(ScatteredField {
                surface,
                incident,
                solution: finish_scattering(surface, solution, p, v, law)?,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source(position_m: [f64; 3]) -> Monopole {
        Monopole {
            position_m,
            volume_velocity_m3_s: C64::new(1e-4, -3e-5),
        }
    }

    fn same_bits(a: &[C64], b: &[C64]) {
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(b) {
            assert_eq!(a.re.to_bits(), b.re.to_bits());
            assert_eq!(a.im.to_bits(), b.im.to_bits());
        }
    }

    /// G0: the exact spherical field obeys Euler's equation, and its static
    /// volume flux has the declared positive outward sign. Omitting 1/(kr)
    /// in the velocity or borrowing the opposite phasor convention fails.
    #[test]
    fn monopole_keeps_reactive_velocity_and_outward_volume_flux() {
        let medium = Medium::air();
        let sources = [source([0.1, -0.2, 0.3])];
        let incident = IncidentField::Monopoles(&sources);
        let point = [0.7, 0.4, 0.6];
        let k = 3.0;
        let sample = incident.pressure_velocity_at(k, medium, point).unwrap();
        for axis in 0..3 {
            let mut plus = point;
            let mut minus = point;
            plus[axis] += 1e-6;
            minus[axis] -= 1e-6;
            let gradient = (incident
                .pressure_velocity_at(k, medium, plus)
                .unwrap()
                .pressure_pa
                - incident
                    .pressure_velocity_at(k, medium, minus)
                    .unwrap()
                    .pressure_pa)
                .scale(0.5e6);
            let momentum =
                sample.velocity_m_s[axis] * C64::new(0.0, k * medium.sound_speed * medium.density);
            assert!((gradient - momentum).abs() < 2e-9 * gradient.abs());
        }
        let radius = 0.02;
        let low = incident
            .pressure_velocity_at(1e-6, medium, [0.1 + radius, -0.2, 0.3])
            .unwrap();
        let flux =
            low.velocity_m_s[0].scale(4.0 * PI * radius * radius) / sources[0].volume_velocity_m3_s;
        assert!((flux - C64::ONE).abs() < 2e-14, "outward flux: {flux:?}");
        let pressure_over_flow = low.pressure_pa / sources[0].volume_velocity_m3_s;
        assert!(pressure_over_flow.im < 0.0);
    }

    #[test]
    fn finite_sources_refuse_interior_surface_and_singular_receiver() {
        let surface = SpherePanels::icosphere(1.0, 0).unwrap();
        let medium = Medium::air();
        for position in [[0.0; 3], surface.centroids()[0]] {
            let sources = [source(position)];
            assert!(
                solve_scattering(
                    &surface,
                    0.6,
                    medium,
                    Admittance::Rigid,
                    IncidentField::Monopoles(&sources),
                    0.01,
                    Formulation::BurtonMiller
                )
                .is_err()
            );
        }
        let sources = [source([0.0, 0.0, 3.0])];
        let field = solve_scattering(
            &surface,
            0.6,
            medium,
            Admittance::Rigid,
            IncidentField::Monopoles(&sources),
            0.01,
            Formulation::BurtonMiller,
        )
        .unwrap();
        assert!(
            field
                .total_pressure_at(
                    &[sources[0].position_m],
                    0.01,
                    near_field::Options::default()
                )
                .is_err()
        );
        let too_many = vec![sources[0]; MAX_MONOPOLES + 1];
        assert!(
            IncidentField::Monopoles(&too_many)
                .pressure_velocity_at(0.6, medium, [0.0, 0.0, 2.0])
                .is_err()
        );
    }

    /// G5: sharing the factorization must not change any boundary result,
    /// including passive Robin work and the historical plane-wave solution.
    #[test]
    fn diffraction_batch_preserves_single_rhs_and_plane_wave_bits() {
        let surface = SpherePanels::icosphere(1.0, 1).unwrap();
        let medium = Medium::air();
        let a = [source([0.0, 0.0, 3.0])];
        let b = [source([2.7, 0.4, 0.0])];
        let plane = PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.7, 0.2),
        };
        let incidents = [
            IncidentField::Monopoles(&a),
            IncidentField::Monopoles(&b),
            IncidentField::PlaneWave(plane),
        ];
        let beta = vec![C64::new(0.4, 0.2); surface.areas().len()];
        for law in [
            Admittance::Rigid,
            Admittance::Local {
                beta: &beta,
                declared_active: false,
            },
        ] {
            let batch = solve_scattering_batch(
                &surface,
                0.7,
                medium,
                law,
                &incidents,
                0.01,
                Formulation::BurtonMiller,
            )
            .unwrap();
            for (result, &incident) in batch.iter().zip(&incidents) {
                let single = solve_scattering(
                    &surface,
                    0.7,
                    medium,
                    law,
                    incident,
                    0.01,
                    Formulation::BurtonMiller,
                )
                .unwrap();
                same_bits(
                    &result.solution.scattered.pressure,
                    &single.solution.scattered.pressure,
                );
                same_bits(
                    &result.solution.scattered.velocity,
                    &single.solution.scattered.velocity,
                );
                assert_eq!(
                    result.solution.absorbed_power_w.to_bits(),
                    single.solution.absorbed_power_w.to_bits()
                );
            }
            let old =
                solve_scattering_bc(&surface, 0.7, medium, law, plane, Formulation::BurtonMiller)
                    .unwrap();
            same_bits(
                &batch[2].solution.scattered.pressure,
                &old.scattered.pressure,
            );
        }
    }

    // Independent spherical-Bessel recurrences for a rigid sphere, not the
    // BEM's Green quadrature. NIST DLMF 10.60.1-2 supply the addition theorem:
    // G=ik/(4pi) sum (2l+1) j_l(kr_<) h_l(kr_>) P_l(cos gamma).
    // https://dlmf.nist.gov/10.60
    // The source is on +z, outside the unit sphere.
    fn spherical(x: f64, maximum: usize) -> (Vec<f64>, Vec<f64>) {
        let top = maximum + 48;
        let mut raw = vec![0.0; top + 2];
        raw[top] = 1.0;
        for n in (1..=top).rev() {
            raw[n - 1] = (2 * n + 1) as f64 / x * raw[n] - raw[n + 1];
        }
        let scale = det::sin(x) / x / raw[0];
        let j: Vec<_> = raw[..=maximum].iter().map(|v| v * scale).collect();
        let mut y = vec![-det::cos(x) / x, -det::cos(x) / (x * x) - det::sin(x) / x];
        for n in 1..maximum {
            y.push((2 * n + 1) as f64 / x * y[n] - y[n - 1]);
        }
        (j, y)
    }

    fn sphere_pressure(k: f64, medium: Medium, source: Monopole, receiver: [f64; 3]) -> C64 {
        let maximum = 24;
        let radius = norm(receiver);
        let mu = receiver[2] / radius;
        let (j, y) = spherical(k, maximum + 1);
        let (js, ys) = spherical(k * source.position_m[2], maximum);
        let (jr, yr) = spherical(k * radius, maximum);
        let mut previous = 1.0;
        let mut legendre = 1.0;
        let mut scattered = C64::ZERO;
        for n in 0..=maximum {
            if n == 1 {
                legendre = mu;
            } else if n > 1 {
                let next =
                    ((2 * n - 1) as f64 * mu * legendre - (n - 1) as f64 * previous) / n as f64;
                previous = legendre;
                legendre = next;
            }
            let jp = if n == 0 {
                -j[1]
            } else {
                j[n - 1] - (n + 1) as f64 / k * j[n]
            };
            let yp = if n == 0 {
                -y[1]
            } else {
                y[n - 1] - (n + 1) as f64 / k * y[n]
            };
            scattered = scattered
                - C64::new(js[n], ys[n]) * C64::new(jr[n], yr[n]) / C64::new(jp, yp)
                    * C64::from_re(jp * (2 * n + 1) as f64 * legendre);
        }
        let incident = IncidentField::Monopoles(&[source])
            .pressure_velocity_at(k, medium, receiver)
            .unwrap()
            .pressure_pa;
        incident
            + scattered
                * source
                    .volume_velocity_m3_s
                    .scale(k * k * medium.sound_speed * medium.density / (4.0 * PI))
    }

    /// G1/G2: an exterior finite source's complex total field converges to the
    /// separated spherical-wave solution on lit, side and shadow receivers.
    #[test]
    fn finite_source_diffraction_converges_to_rigid_sphere_series() {
        let medium = Medium::air();
        let k = 0.7;
        let sources = [source([0.0, 0.0, 3.0])];
        let receivers = [[0.0, 0.0, 2.0], [2.0, 0.0, 0.0], [0.0, 0.0, -2.0]];
        let exact = receivers.map(|x| sphere_pressure(k, medium, sources[0], x));
        let mut errors = Vec::new();
        for rung in 0..=2 {
            let surface = SpherePanels::icosphere(1.0, rung).unwrap();
            let solved = solve_scattering(
                &surface,
                k,
                medium,
                Admittance::Rigid,
                IncidentField::Monopoles(&sources),
                0.01,
                Formulation::PlainCbie,
            )
            .unwrap();
            let pressure = solved
                .total_pressure_at(&receivers, 0.01, near_field::Options::default())
                .unwrap()
                .pressure;
            let error = pressure
                .iter()
                .zip(exact)
                .map(|(&a, b)| (a - b).norm_sq())
                .sum::<f64>()
                .sqrt()
                / exact.iter().map(|p| p.norm_sq()).sum::<f64>().sqrt();
            errors.push(error);
        }
        assert!(
            errors.windows(2).all(|pair| pair[1] < pair[0]),
            "sphere refinement: {errors:?}"
        );
        assert!(
            errors[2] < 0.02,
            "finite-source complex sphere error: {errors:?}"
        );
    }

    /// G3: a reciprocal rigid obstacle cannot depend on which endpoint is
    /// called the emitter. This checks source/receiver signs and diffraction.
    #[test]
    fn finite_source_diffraction_is_reciprocal() {
        let surface = SpherePanels::icosphere(1.0, 2).unwrap();
        let medium = Medium::air();
        let a = [source([0.2, 0.1, 2.3])];
        let b = [source([-2.4, 0.4, -0.7])];
        let incidents = [IncidentField::Monopoles(&a), IncidentField::Monopoles(&b)];
        let fields = solve_scattering_batch(
            &surface,
            0.8,
            medium,
            Admittance::Rigid,
            &incidents,
            0.01,
            Formulation::PlainCbie,
        )
        .unwrap();
        let ab = fields[0]
            .total_pressure_at(&[b[0].position_m], 0.01, near_field::Options::default())
            .unwrap()
            .pressure[0];
        let ba = fields[1]
            .total_pressure_at(&[a[0].position_m], 0.01, near_field::Options::default())
            .unwrap()
            .pressure[0];
        assert!(
            (ab - ba).abs() / ab.abs() < 1e-3,
            "reciprocity: {ab:?}, {ba:?}"
        );
    }

    /// G3: the source and scattered fields share one outgoing phase and one
    /// geometric spreading factor. This also catches a missing direct source.
    #[test]
    fn total_far_field_matches_finite_receiver_limit() {
        let surface = SpherePanels::icosphere(1.0, 1).unwrap();
        let medium = Medium::air();
        let sources = [source([0.0, 0.0, 3.0])];
        let k = 0.7;
        let field = solve_scattering(
            &surface,
            k,
            medium,
            Admittance::Rigid,
            IncidentField::Monopoles(&sources),
            0.01,
            Formulation::PlainCbie,
        )
        .unwrap();
        let direction = [0.0, 0.0, -1.0];
        let radius = 5000.0;
        let expected = field.total_far_field(direction).unwrap();
        let measured = field
            .total_pressure_at(
                &[direction.map(|x| x * radius)],
                0.01,
                near_field::Options::default(),
            )
            .unwrap()
            .pressure[0]
            * expik(-k * radius).scale(radius);
        assert!((measured - expected).abs() / expected.abs() < 1e-3);
    }

    /// G2/G4: changing the solve route preserves the original finite-source
    /// Robin field and reports convergence; a pre-cancelled request does no
    /// source admission or operator work. This small direct-near tree isolates
    /// the source adapter from the separately tested FMM truncation ladder.
    #[test]
    fn iterative_finite_source_matches_dense_and_honors_cancellation() {
        let surface = SpherePanels::icosphere(1.0, 1).unwrap();
        let medium = Medium::air();
        let sources = [source([0.0, 0.0, 3.0])];
        let incident = IncidentField::Monopoles(&sources);
        let beta = vec![C64::new(0.3, -0.1); surface.areas().len()];
        let law = Admittance::Local {
            beta: &beta,
            declared_active: false,
        };
        let options = IterativeOptions {
            tol: 1e-10,
            ..IterativeOptions::default()
        };
        let gate = CancelGate::new_clock_free();
        let dense = solve_scattering(
            &surface,
            0.7,
            medium,
            law,
            incident,
            0.01,
            Formulation::BurtonMiller,
        )
        .unwrap();
        let (iterative, report, receipt) = solve_scattering_iterative(
            &surface,
            0.7,
            medium,
            law,
            incident,
            0.01,
            Formulation::BurtonMiller,
            &options,
            &gate,
        )
        .unwrap();
        assert!(report.converged);
        assert!(receipt.matvecs > 0);
        let exact = &dense.solution.scattered.pressure;
        let error = iterative
            .solution
            .scattered
            .pressure
            .iter()
            .zip(exact)
            .map(|(&a, &b)| (a - b).norm_sq())
            .sum::<f64>()
            .sqrt()
            / exact.iter().map(|p| p.norm_sq()).sum::<f64>().sqrt();
        assert!(error < 1e-8, "finite-source dense/iterative error: {error}");
        gate.request();
        assert!(matches!(
            solve_scattering_iterative(
                &surface,
                f64::NAN,
                medium,
                law,
                incident,
                0.01,
                Formulation::BurtonMiller,
                &options,
                &gate,
            ),
            Err(IterativeError::Cancelled { matvecs: 0 })
        ));
    }
}
