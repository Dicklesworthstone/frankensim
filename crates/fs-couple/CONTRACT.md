

## Reflected pressure at an original moving microphone

`MovingCompactRenderer::with_moving_receiver_reflections` borrows `PreparedReflectedMovingReceiverSource` histories sharing the same original `MovingReflectionObserver`. Every physical strength interpolant is compiled once; each direct/reflected path queries it at its own retarded source time and contributes signed physical pressure to the existing compensated sum. The microphone's trajectory changes arrival times and source-clock bandwidth, never the Eulerian pressure coefficient law.

The existing source/path counts, complete finite-history checks, source/receiver-pair budgets, retained root work, global sample clock and `PressureRenderer` interface remain the owners. Bandwidth admission includes every image's reflected receiver velocity; a direct-path relative-motion bound is insufficient. A foreign microphone/wall owner, incomplete reflected prehistory, a tighter unsupported residual policy or an out-of-window clock refuses before rendering. Callback partitioning remains bit-identical, coherent opposing sources cancel, and cancelled blocks preserve both output and accepted sample count. The underlying stationary ideal-wall, quiescent lossless compact-source model boundaries remain explicit; no trajectory, radiation feedback or PCM encoder semantics are introduced.
