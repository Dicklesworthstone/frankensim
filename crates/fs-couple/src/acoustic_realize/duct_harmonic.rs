//! Periodic radiation from the original duct's complete outlet-flow rows.
//!
//! A radiating mouth need not be an external power port. Its original
//! [`AcousticOutletFlow`] retains loaded-terminal inertia, native orientation,
//! and every material-memory coordinate. Contracting that row against the
//! complete [`HbOrbit`] observes the actual terminal; it does not substitute
//! inlet flow, add another radiation load, or reconstruct a second duct.

use crate::pcm_wav::baffled::{
    BaffledPressure, CircularOutletReceiver, RayleighMedium,
    harmonic::{Field, MAX_KERNEL_TERMS, Options, Prepared},
};
use fs_aeroac::{
    AeroacError,
    observation::{HarmonicFieldLine, HarmonicObservationBudget, PreparedHarmonicField},
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};
use fs_orbit::HbOrbit;
use fs_phs::acoustic_outlet::{AcousticOutletFlow, AcousticTapFlow};

/// Native opening identity; tap indices include closed entries in the source.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DuctOpeningId {
    Bell,
    Tap(usize),
}

/// Borrow the rows returned with the SAME original acoustic-chain assembly.
pub struct DuctOpeningFlows<'a> {
    pub bell: &'a AcousticOutletFlow,
    /// Actual original terminal radius; flow rows do not encode geometry.
    pub bell_radius_m: f64,
    /// Active source taps in their original order, including sparse tap IDs.
    pub taps: &'a [AcousticTapFlow],
}

/// One physical receiver expressed relative to each native opening's centre.
///
/// Positions and normals use one common Cartesian frame. The caller supplies
/// the outward aperture normals; tap stations do not establish hole azimuths
/// or instrument geometry. The original Rayleigh owner receives each opening
/// in its local +z frame; velocity is rotated back before coherent summation.
/// Each opening has its own infinite-baffle image; mutual scattering/loading
/// is not introduced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctOpeningReceiver {
    pub opening: DuctOpeningId,
    /// Finite nonzero outward direction, normalized during admission.
    pub outward_normal: [f64; 3],
    /// Receiver minus opening centre, in the common Cartesian frame.
    pub geometry: CircularOutletReceiver,
}

/// Explicit total work/storage bounds for a cold periodic observation.
#[derive(Clone, Copy, Debug)]
pub struct DuctHarmonicBudget {
    pub max_openings: usize,
    /// All retained state harmonics must fit; none are silently discarded.
    pub max_harmonics: usize,
    /// Maximum openings * original states * (harmonics + 1).
    pub max_projection_terms: usize,
    /// Sum of quadrature points over all openings, times positive harmonics.
    pub max_kernel_terms: usize,
    /// Conservative numeric-array envelope; caller state/rows are borrowed.
    /// Allocator metadata is excluded. Admission precedes numeric allocation.
    pub max_workspace_bytes: usize,
    /// Original Rayleigh receiver's positive separation requirement [m].
    pub minimum_clearance_m: f64,
}

/// One unchanged native mouth's physical spectrum and observation geometry.
#[derive(Clone, Debug)]
pub struct DuctOpeningHarmonics {
    pub opening: DuctOpeningId,
    /// Native axial fraction for a tap; the bell has no tap station.
    pub station: Option<f64>,
    /// Native effective open radius for taps, without a second area multiplier.
    pub radius_m: f64,
    /// Normalized aperture normal in the common Cartesian frame.
    pub outward_normal: [f64; 3],
    pub receiver: CircularOutletReceiver,
    /// Static outward flow [m^3/s], retained but never radiated as an AC line.
    pub dc_volume_flow_m3_s: f64,
    /// Peak volume-flow phasors, index h-1, under exp(-i h omega t).
    pub peak_volume_flow_m3_s: Vec<C64>,
}

/// Coherent receiver field ready for the existing harmonic pressure renderer.
#[derive(Clone, Debug)]
pub struct DuctHarmonicField {
    pub omega_rad_s: f64,
    /// Explicit lossless homogeneous exterior, independent of the interior.
    pub exterior: FluidAcousticState,
    /// Bell first, then active native taps without reordering or renumbering.
    pub openings: Vec<DuctOpeningHarmonics>,
    pub field: PreparedHarmonicField<3>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DuctHarmonicError {
    Invalid(&'static str),
    Budget(&'static str),
    /// Original circular quadrature or harmonic Rayleigh refusal.
    Radiation(String),
    Acoustic(AeroacError),
    Cancelled,
}
impl core::fmt::Display for DuctHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native duct harmonic observation: {self:?}")
    }
}
impl std::error::Error for DuctHarmonicError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DuctHarmonicError> {
    if cancelled() {
        Err(DuctHarmonicError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Observe native loaded bell/tap flow through the existing Rayleigh owner.
///
/// Every original state row participates in the direct DD Fourier contraction,
/// costing O(openings * states * harmonics); no sampled reconstruction or
/// frequency fitting is used. Real DC is retained separately. Positive-time
/// Fourier coefficients become negative-time peak phasors by `2 conj(c_h)`.
/// Source output order remains the state order: these flow maps are linear
/// observations of the complete native state, unlike nonlinear storage ports.
///
/// Each opening reuses `BaffledPressure::circular_surface` at its native radius
/// and `harmonic::Prepared` at every retained positive frequency. Individual
/// fields are combined coherently in DD before the existing field/renderer
/// owner samples them. Absolute source phase and propagation remain intact.
/// The opening receivers must identify one common receiver in a common frame.
/// The caller must supply the rows from the same assembly as the orbit; shape
/// admission cannot establish a source identity. Source objects are unchanged.
///
/// This is a one-way, homogeneous lossless exterior observation with uniform
/// flow over each original circular mouth. It neither validates the candidate
/// orbit, changes its radiation impedance, adds atmospheric absorption, nor
/// resolves mutual aperture loading or scattering between oriented baffles.
/// The original quadrature/band/clearance refusals remain binding. No spatial
/// error enclosure or continuum radiation-match claim is made.
///
/// All work is cold and bounded. Radiation costs O(harmonics * total points),
/// and the existing field preparation is O(harmonics^2). Cancellation is polled
/// during projection and the original receiver preparation; refusal publishes
/// no partial field. Use the existing `HarmonicPressureRenderer` for output.
pub fn observe_duct_harmonics(
    orbit: &HbOrbit,
    flows: DuctOpeningFlows<'_>,
    receivers: &[DuctOpeningReceiver],
    exterior: FluidAcousticState,
    budget: DuctHarmonicBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DuctHarmonicField, DuctHarmonicError> {
    poll(&mut cancelled)?;
    let states = orbit.coeffs.len();
    let width = orbit.coeffs.first().map_or(0, Vec::len);
    let mouths = flows
        .taps
        .len()
        .checked_add(1)
        .ok_or(DuctHarmonicError::Budget("opening count overflow"))?;
    if states == 0
        || width == 0
        || !orbit.omega.is_finite()
        || orbit.omega <= 0.0
        || receivers.len() != mouths
        || !flows.bell_radius_m.is_finite()
        || flows.bell_radius_m <= 0.0
        || !budget.minimum_clearance_m.is_finite()
        || budget.minimum_clearance_m <= 0.0
        || orbit
            .coeffs
            .iter()
            .any(|row| row.len() != width || row[0].im != 0.0 || row.iter().any(|&z| !finite(z)))
    {
        return Err(DuctHarmonicError::Invalid(
            "complete finite state, opening geometry and receiver rows required",
        ));
    }
    let harmonics = width - 1;
    let fundamental = orbit.omega / core::f64::consts::TAU;
    let highest = harmonics as f64 * fundamental;
    if fundamental <= 0.0
        || !highest.is_finite()
        || !((harmonics.max(1) as f64) * orbit.omega).is_finite()
    {
        return Err(DuctHarmonicError::Invalid(
            "physical frequency range is unrepresentable",
        ));
    }
    let projection = mouths
        .checked_mul(states)
        .and_then(|v| v.checked_mul(width));
    let pairs = harmonics
        .checked_mul(harmonics)
        .ok_or(DuctHarmonicError::Budget("field line-pair count overflow"))?;
    if mouths > budget.max_openings
        || harmonics > budget.max_harmonics
        || projection.is_none_or(|n| n > budget.max_projection_terms)
    {
        return Err(DuctHarmonicError::Budget(
            "complete outlet projection exceeds its work budget",
        ));
    }
    let mut points = 0_usize;
    let mut largest_surface = 0_usize;
    for (index, receiver) in receivers.iter().enumerate() {
        poll(&mut cancelled)?;
        let (id, row, radius, _) = opening(&flows, index);
        let (_, geometry) = frame(receiver)?;
        if receiver.opening != id
            || row.len() != states
            || row.iter().any(|v| !v.is_finite())
            || !radius.is_finite()
            || radius <= 0.0
            || geometry.position_m.iter().any(|v| !v.is_finite())
            || geometry.position_m[2] < budget.minimum_clearance_m
            || !(1..=128).contains(&geometry.radial_rings)
            || !(8..=512).contains(&geometry.angular_points)
            || !geometry.maximum_frequency_hz.is_finite()
            || geometry.maximum_frequency_hz <= 0.0
            || highest > geometry.maximum_frequency_hz
        {
            return Err(DuctHarmonicError::Invalid(
                "native opening identity, complete flow row or resolved receiver band disagrees",
            ));
        }
        if index > 0 {
            let tap = &flows.taps[index - 1];
            if !tap.station.is_finite()
                || !(0.0..=1.0).contains(&tap.station)
                || flows.taps[..index - 1]
                    .iter()
                    .any(|other| other.tap_index == tap.tap_index)
            {
                return Err(DuctHarmonicError::Invalid(
                    "native tap IDs must be distinct and stations valid",
                ));
            }
        }
        let count = geometry.radial_rings * geometry.angular_points;
        points = points
            .checked_add(count)
            .ok_or(DuctHarmonicError::Budget("quadrature count overflow"))?;
        largest_surface = largest_surface.max(count);
    }
    if points
        .checked_mul(harmonics)
        .is_none_or(|n| n > budget.max_kernel_terms)
    {
        return Err(DuctHarmonicError::Budget(
            "all opening-frequency kernels must fit the work budget",
        ));
    }
    // Flow spectra, current quadrature/shapes, field lines, DD accumulators,
    // coalesced fields and the original receiver's bounded stack scratch.
    let bytes = 16 * mouths as u128 * harmonics as u128
        + 256 * mouths as u128
        + 512 * width as u128
        + 128 * largest_surface as u128
        + 65_536;
    if bytes > budget.max_workspace_bytes as u128 {
        return Err(DuctHarmonicError::Budget(
            "native outlet observation exceeds numeric storage budget",
        ));
    }
    let medium = RayleighMedium {
        density: exterior.density_kg_m3(),
        sound_speed: exterior.sound_speed_m_s(),
    };
    let mut sums = vec![[[Dd::ZERO; 2]; 4]; harmonics];
    let mut observed = Vec::with_capacity(mouths);
    for (index, receiver) in receivers.iter().enumerate() {
        poll(&mut cancelled)?;
        let (id, row, radius, station) = opening(&flows, index);
        let (basis, geometry) = frame(receiver)?;
        let mut spectrum = Vec::with_capacity(width);
        for harmonic in 0..width {
            let mut value = [Dd::ZERO; 2];
            for (coordinate, (&weight, state)) in row.iter().zip(&orbit.coeffs).enumerate() {
                if coordinate % 64 == 0 {
                    poll(&mut cancelled)?;
                }
                let weight = Dd::from_f64(weight);
                value[0] = value[0] + weight * Dd::from_f64(state[harmonic].re);
                value[1] = value[1] + weight * Dd::from_f64(state[harmonic].im);
            }
            if harmonic != 0 {
                value[0] = value[0] * Dd::from_f64(2.0);
                value[1] = value[1] * Dd::from_f64(-2.0);
            }
            spectrum.push(complex(value)?);
        }
        let surface = BaffledPressure::circular_surface(radius, geometry, medium)
            .map_err(DuctHarmonicError::Radiation)?;
        for (h, &flow) in spectrum[1..].iter().enumerate() {
            poll(&mut cancelled)?;
            let mut stopped = false;
            let transfer = Prepared::new(
                &surface,
                1,
                &[geometry.position_m],
                (h + 1) as f64 * fundamental,
                medium,
                Options {
                    minimum_clearance_m: budget.minimum_clearance_m,
                    maximum_kernel_terms: budget.max_kernel_terms.min(MAX_KERNEL_TERMS),
                },
                || {
                    let stop = cancelled();
                    stopped |= stop;
                    stop
                },
            )
            .map_err(|error| {
                if stopped {
                    DuctHarmonicError::Cancelled
                } else {
                    DuctHarmonicError::Radiation(error)
                }
            })?;
            let mut field = [Field::default()];
            transfer
                .evaluate_into(&[flow], &mut field)
                .map_err(DuctHarmonicError::Radiation)?;
            let local = field[0].particle_velocity_m_s;
            let mut values = [C64::ZERO; 4];
            values[0] = field[0].pressure_pa;
            for component in 0..3 {
                let mut velocity = [Dd::ZERO; 2];
                for axis in 0..3 {
                    let rotation = Dd::from_f64(basis[axis][component]);
                    velocity[0] = velocity[0] + rotation * Dd::from_f64(local[axis].re);
                    velocity[1] = velocity[1] + rotation * Dd::from_f64(local[axis].im);
                }
                values[component + 1] = complex(velocity)?;
            }
            for (sum, value) in sums[h].iter_mut().zip(values) {
                sum[0] = sum[0] + Dd::from_f64(value.re);
                sum[1] = sum[1] + Dd::from_f64(value.im);
            }
        }
        let dc = spectrum.remove(0).re;
        observed.push(DuctOpeningHarmonics {
            opening: id,
            station,
            radius_m: radius,
            outward_normal: basis[2],
            receiver: receiver.geometry,
            dc_volume_flow_m3_s: dc,
            peak_volume_flow_m3_s: spectrum,
        });
    }
    let mut lines = Vec::with_capacity(harmonics);
    for (h, sum) in sums.into_iter().enumerate() {
        poll(&mut cancelled)?;
        lines.push(HarmonicFieldLine {
            frequency_hz: (h + 1) as f64 * fundamental,
            pressure_pa: complex(sum[0])?,
            particle_velocity_m_s: [complex(sum[1])?, complex(sum[2])?, complex(sum[3])?],
        });
    }
    let field = PreparedHarmonicField::new_cancellable(
        &lines,
        HarmonicObservationBudget {
            max_lines: harmonics,
            max_line_pairs: pairs,
        },
        &mut cancelled,
    )
    .map_err(|error| match error {
        AeroacError::Cancelled { .. } => DuctHarmonicError::Cancelled,
        other => DuctHarmonicError::Acoustic(other),
    })?;
    Ok(DuctHarmonicField {
        omega_rad_s: orbit.omega,
        exterior,
        openings: observed,
        field,
    })
}

// Deterministic orthonormal axes, with the original x/y axes for normal +z.
// Scaling before normalization avoids overflow or underflow for finite normals.
fn frame(
    receiver: &DuctOpeningReceiver,
) -> Result<([[f64; 3]; 3], CircularOutletReceiver), DuctHarmonicError> {
    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }
    fn unit(v: [f64; 3]) -> Option<[f64; 3]> {
        if v.iter().any(|x| !x.is_finite()) {
            return None;
        }
        let scale = v.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        if scale == 0.0 {
            return None;
        }
        let scaled = v.map(|x| x / scale);
        let length = det::sqrt(scaled.iter().map(|x| x * x).sum());
        Some(scaled.map(|x| x / length))
    }
    let normal = unit(receiver.outward_normal).ok_or(DuctHarmonicError::Invalid(
        "finite nonzero outward aperture normal required",
    ))?;
    let axis = if normal[1].abs() <= normal[0].abs() && normal[1].abs() <= normal[2].abs() {
        [0.0, 1.0, 0.0]
    } else if normal[0].abs() <= normal[2].abs() {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let tangent =
        unit(cross(axis, normal)).expect("least-aligned unit axis has a nonzero cross product");
    let basis = [tangent, cross(normal, tangent), normal];
    let mut geometry = receiver.geometry;
    geometry.position_m = basis.map(|axis| {
        axis.iter()
            .zip(receiver.geometry.position_m)
            .map(|(a, x)| a * x)
            .sum()
    });
    Ok((basis, geometry))
}

fn opening<'a>(
    flows: &DuctOpeningFlows<'a>,
    index: usize,
) -> (DuctOpeningId, &'a [f64], f64, Option<f64>) {
    if index == 0 {
        (
            DuctOpeningId::Bell,
            flows.bell.coefficients(),
            flows.bell_radius_m,
            None,
        )
    } else {
        let tap = &flows.taps[index - 1];
        (
            DuctOpeningId::Tap(tap.tap_index),
            tap.flow.coefficients(),
            tap.effective_radius_m,
            Some(tap.station),
        )
    }
}

fn complex(value: [Dd; 2]) -> Result<C64, DuctHarmonicError> {
    if value.iter().any(|v| !v.hi.is_finite() || !v.lo.is_finite()) {
        return Err(DuctHarmonicError::Invalid(
            "physical outlet flow or coherent field is unrepresentable",
        ));
    }
    let result = C64::new(value[0].to_f64(), value[1].to_f64());
    if !finite(result) {
        return Err(DuctHarmonicError::Invalid(
            "physical outlet flow or coherent field is unrepresentable",
        ));
    }
    Ok(result)
}
