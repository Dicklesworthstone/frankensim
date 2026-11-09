//! Complete finite-time reflection between two infinite, opposing ideal walls.
//!
//! In local normal coordinates 0 < y < L the source images are y+2kL and
//! -y+2kL, k in Z. Their pressure characters are (s0*sL)^k and
//! s0*(s0*sL)^k. Each reflected force and stress retains its complete vector
//! or tensor transform through the ORIGINAL source kernel at a virtual
//! receiver; its velocity then returns to the physical frame. This is the
//! ideal-boundary image construction of Allen & Berkley, JASA 65 (1979),
//! 943–950, <https://doi.org/10.1121/1.382599>.
//!
//! A causal observation of duration T needs only finitely many images:
//! an omitted |k| >= N+1 has path length >= 2*N*L, since source and receiver
//! both lie between the walls. Preparation selects N > c*T/(2*L), with a
//! guard period for rounding. The omitted arrivals are therefore outside the
//! declared window, not a reflection-order approximation. Initial exterior
//! and source fields are zero. Queries beyond that window refuse.
//!
//! This is a stationary, homogeneous, lossless fluid with exact rigid or
//! pressure-release walls, not a constant-coefficient model of absorbing
//! walls. Finite impedance, wall motion, finite edges, general room geometry,
//! molecular absorption, source feedback and experimental qualification are
//! outside this owner. Source compactness and its history remain caller-owned.

use super::{CompactTransfer, RetardedSourceField, initially_resting, invalid, poll};
use super::super::super::Transfer;
use crate::AeroacError;
use crate::curle3d::array::wedge::reflection;
use crate::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::dd::Dd;

pub use super::TransientAcousticSource;

const SCOPE: &str = "Initially resting prescribed stationary compact 3D volume, force and \
stress sources in a homogeneous lossless fluid between two infinite parallel ideal rigid \
or pressure-release walls. Every image arriving within the explicit finite window retains \
its own delay and full coherent pressure/vector velocity. No finite impedance, finite edges, \
motion, absorption, feedback, source/bandwidth validation or measured-SPL qualification.";

/// Two opposing parallel walls in one original Cartesian frame.
/// The lower plane's inward normal points toward the upper wall a distance L
/// away. The separation is authored directly; no supplied angle is snapped.
#[derive(Clone, Copy, Debug)]
pub struct IdealAcousticSlab {
    lower: IdealAcousticPlane,
    separation_m: f64,
    upper_boundary: IdealPlaneBoundary,
}

impl IdealAcousticSlab {
    pub fn new(
        lower: IdealAcousticPlane,
        separation_m: f64,
        upper_boundary: IdealPlaneBoundary,
    ) -> Result<Self, AeroacError> {
        if !separation_m.is_finite() || separation_m <= 0.0 {
            return Err(invalid("acoustic slab separation must be positive and finite"));
        }
        Ok(Self { lower, separation_m, upper_boundary })
    }

    pub const fn lower_plane(&self) -> IdealAcousticPlane { self.lower }
    pub const fn separation_m(&self) -> f64 { self.separation_m }
    pub const fn upper_boundary(&self) -> IdealPlaneBoundary { self.upper_boundary }

    fn admit_position(&self, position: [f64; 3], source: bool) -> Result<(), AeroacError> {
        if position.iter().any(|v| !v.is_finite()) {
            return Err(invalid("acoustic slab positions must be finite"));
        }
        let distance = self.lower.distance(position)?;
        let upper = Dd::from_f64(self.separation_m);
        let inside = if source {
            Dd::ZERO.lt(distance) && distance.lt(upper)
        } else {
            !distance.lt(Dd::ZERO) && !upper.lt(distance)
        };
        if !inside {
            return Err(invalid("source must be inside the slab; receiver may be on its walls"));
        }
        Ok(())
    }

    // Invert the source image on the receiver, keeping the original source
    // moments and physical clock. This avoids per-query vector/tensor images.
    fn virtual_receiver(
        &self,
        receiver: [f64; 3],
        period: i64,
        mirrored: bool,
    ) -> Result<[f64; 3], AeroacError> {
        let shift = Dd::from_f64(2.0 * period as f64) * Dd::from_f64(self.separation_m);
        let offset = if mirrored {
            shift - self.lower.distance(receiver)? * Dd::from_f64(2.0)
        } else {
            -shift
        };
        let normal = self.lower.inward_normal();
        let image = core::array::from_fn(|axis| {
            (Dd::from_f64(receiver[axis]) + offset * Dd::from_f64(normal[axis])).to_f64()
        });
        if image.iter().any(|v| !v.is_finite()) {
            return Err(invalid("acoustic slab virtual receiver is unrepresentable"));
        }
        Ok(image)
    }
}

/// A causal source epoch and a finite observation duration [s].
/// Every source and the exterior are initially at rest before the epoch.
/// A separate duration retains short windows even at a large absolute epoch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CausalSlabWindow {
    pub source_epoch_s: f64,
    pub duration_s: f64,
}

struct SourcePath {
    source: usize,
    field: Transfer,
    compact: CompactTransfer,
}

/// Original sources and every image needed for the declared finite window.
/// Geometry, source parity and all near/far coefficients are cached once.
/// Evaluation has O(paths) work, no allocation or transcendentals, and no
/// source callback before a path arrives. There is no reflection-order knob.
pub struct PreparedSlabSourceObserver {
    slab: IdealAcousticSlab,
    acoustic: FluidAcousticState,
    receiver_m: [f64; 3],
    window: CausalSlabWindow,
    source_count: usize,
    paths: Vec<SourcePath>,
}

impl PreparedSlabSourceObserver {
    /// Prepare the complete finite-time image lattice, including direct paths.
    /// `maximum_source_paths` bounds ALL retained image/source pairs, not just
    /// audible paths or physical sources. Insufficient work or allocation
    /// budgets refuse instead of dropping later echoes. Empty arrays are silence.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        slab: IdealAcousticSlab,
        acoustic: FluidAcousticState,
        receiver_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        window: CausalSlabWindow,
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(&mut cancelled)?;
        if !window.source_epoch_s.is_finite()
            || !window.duration_s.is_finite()
            || window.duration_s <= 0.0
        {
            return Err(invalid("acoustic slab needs a finite epoch and positive duration"));
        }
        slab.admit_position(receiver_m, false)?;
        for &source in source_positions_m {
            poll(&mut cancelled)?;
            slab.admit_position(source, true)?;
        }
        let one_way_s = slab.separation_m / acoustic.sound_speed_m_s();
        let periods = (window.duration_s / one_way_s) * 0.5;
        // Below 2^48 the quotient's accumulated f64 roundoff is smaller than
        // the extra full period. It also keeps every signed image index exact.
        if !one_way_s.is_finite() || one_way_s <= 0.0 || !periods.is_finite()
            || periods > (1_u64 << 48) as f64
        {
            return Err(invalid("acoustic slab propagation extent is unrepresentable"));
        }
        let n = (periods.ceil() as usize).checked_add(1)
            .ok_or_else(|| invalid("acoustic slab image index extent overflows"))?;
        let count = n.checked_mul(4).and_then(|n| n.checked_add(2))
            .and_then(|n| n.checked_mul(source_positions_m.len()))
            .filter(|&n| n <= maximum_source_paths)
            .ok_or_else(|| invalid("complete acoustic slab paths exceed preparation budget"))?;
        let mut paths = Vec::new();
        paths.try_reserve_exact(count)
            .map_err(|_| invalid("acoustic slab paths exceed available memory"))?;
        let mut result = Self {
            slab, acoustic, receiver_m, window,
            source_count: source_positions_m.len(), paths,
        };
        // Preserve original direct-source order before any reflected group.
        result.append_image(0, false, source_positions_m, &mut cancelled)?;
        if result.source_count != 0 {
            for period in -(n as i64)..=n as i64 {
                poll(&mut cancelled)?;
                if period != 0 {
                    result.append_image(period, false, source_positions_m, &mut cancelled)?;
                }
                result.append_image(period, true, source_positions_m, &mut cancelled)?;
            }
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn append_image(
        &mut self,
        period: i64,
        mirrored: bool,
        sources: &[[f64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AeroacError> {
        let receiver = self.slab.virtual_receiver(self.receiver_m, period, mirrored)?;
        let mixed = self.slab.lower.boundary() != self.slab.upper_boundary;
        let mut sign = if mixed && period % 2 != 0 { -1.0 } else { 1.0 };
        if mirrored { sign *= self.slab.lower.boundary().sign(); }
        let matrix = mirrored.then(|| reflection(self.slab.lower.inward_normal()));
        for (source, &position) in sources.iter().enumerate() {
            poll(cancelled)?;
            let mut field = Transfer::new(self.acoustic, receiver, position)?;
            field.pressure_scales = field.pressure_scales.map(|v| v * sign);
            field.velocity_scales = field.velocity_scales.map(|v| v * sign);
            field.inverse_image = matrix;
            let compact = CompactTransfer::new(&field, self.acoustic)?;
            self.paths.push(SourcePath { source, field, compact });
        }
        Ok(())
    }

    pub const fn slab(&self) -> IdealAcousticSlab { self.slab }
    pub const fn acoustic(&self) -> FluidAcousticState { self.acoustic }
    pub const fn receiver_position_m(&self) -> [f64; 3] { self.receiver_m }
    pub const fn window(&self) -> CausalSlabWindow { self.window }
    pub const fn source_count(&self) -> usize { self.source_count }
    pub fn path_count(&self) -> usize { self.paths.len() }
    pub fn path_source_index(&self, path: usize) -> Option<usize> {
        self.paths.get(path).map(|p| p.source)
    }
    pub fn path_delay_s(&self, path: usize) -> Option<f64> {
        self.paths.get(path).map(|p| p.field.delay_s)
    }
    pub fn source_position_m(&self, source: usize) -> Option<[f64; 3]> {
        (source < self.source_count).then(|| self.paths[source].field.source_m)
    }

    /// Observe one time `source_epoch + elapsed_s`, within the prepared window.
    ///
    /// Calling declares ZERO prehistory, including the exterior near field.
    /// Initial Q, F, force impulse, stress, stress rate and stress integral
    /// are checked at the epoch for every physical source. Finite rate and
    /// acceleration jumps are allowed, as in the original causal source owner.
    /// This is not permission to silently zero-pad a noncausal record.
    ///
    /// The callback receives the ORIGINAL source index and DD emission time;
    /// it supplies consistent values, derivatives and integrals. Images never
    /// advance or clone its state. Callbacks may occur out of time order and
    /// should be read-only; callback side effects cannot be rolled back.
    /// Pressure and full velocity sum before p*u, retaining coherent echoes
    /// and reactive near fields. Refusal/cancellation publishes no field.
    pub fn evaluate_initially_resting(
        &self,
        elapsed_s: f64,
        mut source_at: impl FnMut(usize, Dd) -> Result<TransientAcousticSource, AeroacError>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RetardedSourceField, AeroacError> {
        poll(&mut cancelled)?;
        if !elapsed_s.is_finite() || elapsed_s < 0.0 || elapsed_s > self.window.duration_s {
            return Err(invalid("acoustic slab query exceeds its complete causal window"));
        }
        let epoch = Dd::from_f64(self.window.source_epoch_s);
        for index in 0..self.source_count {
            poll(&mut cancelled)?;
            if !initially_resting(&source_at(index, epoch)?) {
                return Err(invalid("acoustic slab requires compatible initially resting sources"));
            }
        }
        let mut sum = [Dd::ZERO; 4];
        for (index, path) in self.paths.iter().enumerate() {
            if index.is_multiple_of(256) { poll(&mut cancelled)?; }
            if elapsed_s < path.field.delay_s { continue; }
            let relative = Dd::from_f64(elapsed_s) - Dd::from_f64(path.field.delay_s);
            let emission = epoch + relative;
            if !emission.hi.is_finite() || !emission.lo.is_finite()
                || (relative != Dd::ZERO && !epoch.lt(emission))
            {
                return Err(invalid("acoustic slab emission time is unresolved on the DD clock"));
            }
            let rows = path.compact.inputs(&path.field, source_at(path.source, emission)?)?;
            for row in rows {
                for channel in 0..4 {
                    sum[channel] = sum[channel] + Dd::from_f64(row[channel]);
                }
            }
        }
        poll(&mut cancelled)?;
        let [pressure_pa, ux, uy, uz] = sum.map(Dd::to_f64);
        let particle_velocity_m_s = [ux, uy, uz];
        let instantaneous_intensity_w_m2 = particle_velocity_m_s.map(|u| pressure_pa * u);
        if !pressure_pa.is_finite() || particle_velocity_m_s.iter()
            .chain(&instantaneous_intensity_w_m2).any(|v| !v.is_finite())
        {
            return Err(invalid("coherent acoustic slab field is unrepresentable"));
        }
        Ok(RetardedSourceField {
            pressure_pa, particle_velocity_m_s, instantaneous_intensity_w_m2, scope: SCOPE,
        })
    }
}
