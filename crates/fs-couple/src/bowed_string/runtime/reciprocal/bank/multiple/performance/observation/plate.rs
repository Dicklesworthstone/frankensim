//! Finite-distance observation of the original regional material plate.
//!
//! The same triangle quadrature and signed mass-normalized velocity traces
//! supply the mechanical Rayleigh load and this microphone. The shared piano
//! and wind receiver owns propagation, material-gas loss and dispersion; no
//! compact-area reduction, second force solve or fitted output gain occurs.
use super::{
    BoundReceiver, BoundReport, BoundShellObserver, MultiBowPerformance, ShellObserverError,
    invalid, poll,
};
use crate::baffled_pressure::gas::{GasPropagation, GasPropagationLimits, GasPropagationReport};
use crate::baffled_pressure::{BaffledPressure, SurfaceSample};
use crate::bowed_string::runtime::reciprocal::bank::radiation::rayleigh::surface_samples;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;

/// One stationary receiver in the positive half-space of the original flat
/// plate (z=0). Its numerical histories begin at rest on the mechanical clock.
#[derive(Clone, Debug)]
pub struct PlateObserverOptions {
    pub position_m: [f64; 3],
    /// Declared resolved band for the original surface-resolution diagnostic.
    /// This is not a filter or proof that nonlinear motion stays in that band.
    pub maximum_frequency_hz: f64,
    pub maximum_surface_cells: usize,
    /// Complete retained convolution and derivative multiplies per sample.
    pub maximum_replay_terms: usize,
    pub propagation: GasPropagation,
}
impl PlateObserverOptions {
    /// Cheap complete declaration admission, before source or spectral work.
    pub fn validate(&self, rate: u32) -> Result<(), ShellObserverError> {
        let ceiling = GasPropagationLimits::surface_reference();
        let limits = self.propagation.limits;
        if !(8_000..=192_000).contains(&rate)
            || self.position_m.iter().any(|v| !v.is_finite())
            || self.position_m[2] < 0.05
            || !self.maximum_frequency_hz.is_finite()
            || self.maximum_frequency_hz <= 0.0
            || self.maximum_frequency_hz > 0.45 * f64::from(rate)
            || !(3..=4096).contains(&self.maximum_surface_cells)
            || !(1..=fs_vfit::observation::MAX_OBSERVATION_ENTRIES)
                .contains(&self.maximum_replay_terms)
            || !(1..=ceiling.max_unique_paths).contains(&limits.max_unique_paths)
            || !(1..=ceiling.max_filter_taps).contains(&limits.max_filter_taps)
            || !(1..=ceiling.max_kernel_terms).contains(&limits.max_kernel_terms)
        {
            return Err(invalid(
                "plate observer needs a resolved clock/band, receiver >=5 cm above z=0, and bounded surface/propagation/replay work",
            ));
        }
        Ok(())
    }
}

/// Physical geometry and gas of the actual receiver, including all numerical
/// delay. The load and observer retain their own numerical error limitations.
#[derive(Clone, Debug)]
pub struct PlateObserverReport {
    pub options: PlateObserverOptions,
    pub specimen_identity: ContentHash,
    pub basis_identity: ContentHash,
    pub body_modes: usize,
    pub surface_cells: usize,
    pub minimum_equivalent_disc_ppw: f64,
    pub replay_terms: usize,
    /// FIR support, including physical travel and interpolation taps.
    pub delay_samples: (usize, usize),
    /// Numerical derivative/interpolation/padding delay, excluding flight time.
    pub filter_delay_samples: f64,
    pub gas: GasPropagationReport,
}

impl MultiBowPerformance {
    /// Replace compact pressure with the original plate's finite-distance
    /// Rayleigh microphone and material gas propagation on EACH signed path.
    /// The computed mechanical air reaction remains attached and unchanged.
    ///
    /// This observes one exterior half-space, even when the load reacts from
    /// two independent half-spaces: a second load face is not a second audible
    /// copy. Cavity mouths require their own complete spatial observation and
    /// are refused here rather than silently dropping their sound. Closed
    /// shells, raw fitted loads and an already-started receiver also refuse.
    /// Molecular laws must carry their evaluated source state. Native records
    /// bind the shared molecular declaration to the original performance gas;
    /// an explicit classical selection makes no molecular-loss claim.
    pub fn with_material_plate_observer(
        mut self,
        options: PlateObserverOptions,
        gate: &CancelGate,
    ) -> Result<Self, ShellObserverError> {
        poll(gate)?;
        if self.poisoned
            || self.shell_observer.is_some()
            || self.base_sample != 0
            || self.bank.system.samples_rendered() != 0
        {
            return Err(invalid(
                "plate observation requires one unstarted physical history",
            ));
        }
        let source = self.bank.source.bank();
        let rate = source.options().sample_rate_hz;
        options.validate(rate)?;
        if source.material_cavity().is_some() || source.spatial_cavity().is_some() {
            return Err(invalid(
                "plate observation needs a complete exterior without unobserved cavity mouths",
            ));
        }
        let binding = source.body_radiation().ok_or_else(|| {
            invalid("plate observation requires its original computed Rayleigh air load")
        })?;
        let loading = binding.material_plate_rayleigh().ok_or_else(|| {
            invalid("raw or shell loads have no original plate observation binding")
        })?;
        let body = source.material_chart.as_ref().ok_or_else(|| {
            invalid("plate observation requires its original regional body chart")
        })?;
        let basis = source.material_chart_basis.as_ref().ok_or_else(|| {
            invalid("plate observation requires its original signed body eigenbasis")
        })?;
        let body_component = source.strings().len();
        let modes = self
            .bank
            .system
            .components()
            .get(body_component)
            .ok_or_else(|| invalid("plate observation lost its original body component"))?
            .modes();
        if basis.source_identity != body.specimen.specimen_identity()
            || loading.specimen_identity != basis.source_identity
            || basis.mass_normalized_modes.len() != modes.len()
            || loading.body_modes != modes.len()
            || binding.gas() != source.gas()
            || &loading.gas != source.gas()
            || !binding
                .source_body_coordinates()
                .iter()
                .copied()
                .eq(0..modes.len())
            || modes.iter().any(|m| {
                m.angular_frequency_rad_s > core::f64::consts::TAU * options.maximum_frequency_hz
            })
        {
            return Err(invalid(
                "plate observation must retain every original material mode, gas, load coordinate and resolved reference frequency",
            ));
        }
        let cells = body
            .specimen
            .chart()
            .mesh
            .tris
            .len()
            .checked_mul(3)
            .ok_or_else(|| invalid("plate observer surface extent overflow"))?;
        if cells != loading.surface_cells || cells > options.maximum_surface_cells {
            return Err(invalid(
                "plate observer must retain the complete loaded surface within its cell budget",
            ));
        }
        let samples = surface_samples(&body.specimen.chart().mesh, basis, gate)
            .map_err(|e| ShellObserverError::Plate(e.to_string()))?;
        // Identical surface and gas: the original diagnostic is proportional
        // to wavelength. Reuse it at this declared band instead of allocating
        // another quadratic table of cell-pair distances just for a scalar.
        let minimum_equivalent_disc_ppw = loading.minimum_equivalent_disc_ppw
            * (loading.options.maximum_hz / options.maximum_frequency_hz);
        if !minimum_equivalent_disc_ppw.is_finite() || minimum_equivalent_disc_ppw <= 0.0 {
            return Err(invalid(
                "plate observation wavelength exceeds finite numerical range",
            ));
        }
        let samples: Vec<_> = samples
            .into_iter()
            .map(|p| SurfaceSample {
                position_m: [p.position_m[0], p.position_m[1], 0.0],
                area_m2: p.area_m2,
                mode_shape: p.mode_shape,
            })
            .collect();
        let (mut observer, gas) = BaffledPressure::new_multirate_gas_differentiating(
            &samples,
            modes.len(),
            rate,
            rate,
            options.position_m,
            source.gas(),
            &options.propagation,
            &mut || gate.is_requested(),
        )
        .map_err(ShellObserverError::Plate)?;
        // One clock and zero prehistory let the SAME shared derivative commute
        // through the linear sum, without repeating it for every body mode.
        observer
            .differentiate_pressure()
            .map_err(ShellObserverError::Plate)?;
        let replay_terms = observer.multiply_adds_per_sample();
        if replay_terms > options.maximum_replay_terms {
            return Err(invalid(
                "complete plate observation exceeds its declared replay budget",
            ));
        }
        poll(gate)?;
        let report = PlateObserverReport {
            options,
            specimen_identity: basis.source_identity,
            basis_identity: binding.basis_identity(),
            body_modes: modes.len(),
            surface_cells: cells,
            minimum_equivalent_disc_ppw,
            replay_terms,
            delay_samples: observer.delay_samples,
            filter_delay_samples: observer.filter_delay_samples(),
            gas,
        };
        self.shell_observer = Some(BoundShellObserver {
            report: BoundReport::Plate(Box::new(report)),
            observer: BoundReceiver::Plate(observer),
            coordinates: (0..modes.len())
                .map(|mode| (body_component, mode, 1.0))
                .collect(),
            velocity: vec![0.0; modes.len()],
            maximum_abs_pressure_pa: source.options().network.maximum_abs_pressure_pa,
            samples: 0,
        });
        Ok(self)
    }

    /// Original receiver declaration, signed surface and actual gas-law paths.
    pub fn material_plate_observer(&self) -> Option<&PlateObserverReport> {
        match self
            .shell_observer
            .as_ref()
            .map(|observer| &observer.report)
        {
            Some(BoundReport::Plate(report)) => Some(report),
            _ => None,
        }
    }
}
