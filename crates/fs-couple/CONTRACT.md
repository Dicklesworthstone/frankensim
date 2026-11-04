

## Reflected moving-source histories on the original audio clock

`MovingCompactRenderer::with_ideal_reflections` consumes the original borrowed `PreparedReflectedSourceTrajectory` owners and accepted Q/F histories. Every source shares one original `IdealReflectionObserver`; a separately bound microphone or boundary owner refuses. The existing cubic strength interpolant is compiled once per physical source, then queried at each direct/reflected path's own retarded time. Pressures add coherently in the existing compensated sum before decimation or PCM conversion. No source history or trajectory is copied to create images.

`source_count` continues to count physical sources; `source_path_count` reports all propagation paths. The existing `maximum_source_receiver_pairs` budget counts every path/output pair, and retained root-iteration bounds include that multiplicity. Full history coverage and Doppler bandwidth admission apply separately to every path. Missing reflected prehistory or tails refuse the complete observation instead of becoming silence. Existing global sample indices, no-allocation callbacks, cancellation rollback, pressure units and `PressureRenderer` sinks remain the owners.

Focused G1/G3 tests reconstruct the reflected waveform from independent source-image polynomials, verify original physical-source accounting and exact callback partitioning, retain coherent source cancellation, and exercise path-work, history, boundary identity and transactional cancellation refusals. The underlying ideal stationary-wall, quiescent lossless compact-source model limits remain; this adds no source feedback or experimental acoustic qualification.
