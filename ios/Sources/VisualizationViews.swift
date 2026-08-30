import SwiftUI

struct NativeSimulationCanvas: View {
    let result: SimulationResult?
    let accent: Color
    let isRunning: Bool
    @Environment(\.accessibilityReduceMotion) private var reduceMotion

    var body: some View {
        TimelineView(.animation(minimumInterval: reduceMotion ? 1 : 1 / 24, paused: result == nil || reduceMotion)) { timeline in
            Canvas(opaque: true, colorMode: .linear, rendersAsynchronously: true) { context, size in
                context.fill(Path(CGRect(origin: .zero, size: size)), with: .color(ForgeTheme.background))
                drawGrid(in: &context, size: size)
                if let result {
                    let frame = displayedFrame(for: result, date: timeline.date)
                    switch result.shape {
                    case .grid, .gridFrames:
                        drawField(result, frame: frame, context: &context, size: size)
                    case .xyzPath:
                        drawSpatialPath(result, date: timeline.date, context: &context, size: size)
                    case .triangles:
                        drawTriangles(result, date: timeline.date, context: &context, size: size)
                    case .signal:
                        drawSignal(result, context: &context, size: size)
                    case .campaign:
                        drawCampaign(result, context: &context, size: size)
                    }
                } else if !isRunning {
                    drawIdle(in: &context, size: size)
                }
            }
        }
        .overlay(alignment: .topLeading) {
            HStack(spacing: 7) {
                Circle().fill(isRunning ? ForgeTheme.amber : accent).frame(width: 7, height: 7)
                Text(isRunning ? "KERNEL RUNNING" : result == nil ? "READY" : "LIVE RESULT")
            }
            .font(.system(size: ForgeTheme.size(10), weight: .bold, design: .monospaced))
            .foregroundStyle(isRunning ? ForgeTheme.amber : accent)
            .padding(.horizontal, 10)
            .padding(.vertical, 8)
            .background(.black.opacity(0.54), in: Capsule())
            .padding(12)
        }
        .overlay {
            if isRunning {
                KernelActivityOverlay(accent: accent)
                    .allowsHitTesting(false)
            }
        }
        .clipShape(RoundedRectangle(cornerRadius: 24, style: .continuous))
        .overlay(RoundedRectangle(cornerRadius: 24).stroke(accent.opacity(0.42), lineWidth: 1))
        .shadow(color: accent.opacity(0.13), radius: 28)
        .accessibilityElement(children: .ignore)
        .accessibilityLabel(result == nil ? "Simulation canvas ready" : "Rendered native simulation result")
    }

    private func displayedFrame(for result: SimulationResult, date: Date) -> Int {
        guard result.frames > 1, !reduceMotion else { return max(0, result.frames - 1) }
        return Int(date.timeIntervalSinceReferenceDate * 9) % result.frames
    }

    private func drawGrid(in context: inout GraphicsContext, size: CGSize) {
        var grid = Path()
        let spacing: CGFloat = 36
        stride(from: CGFloat.zero, through: size.width, by: spacing).forEach { x in
            grid.move(to: CGPoint(x: x, y: 0)); grid.addLine(to: CGPoint(x: x, y: size.height))
        }
        stride(from: CGFloat.zero, through: size.height, by: spacing).forEach { y in
            grid.move(to: CGPoint(x: 0, y: y)); grid.addLine(to: CGPoint(x: size.width, y: y))
        }
        context.stroke(grid, with: .color(accent.opacity(0.045)), lineWidth: 0.6)
    }

    private func drawField(_ result: SimulationResult, frame: Int, context: inout GraphicsContext, size: CGSize) {
        guard result.width > 0, result.height > 0 else { return }
        let fieldSize = result.width * result.height
        let start = min(result.values.count, frame * fieldSize)
        let end = min(result.values.count, start + fieldSize)
        guard end - start == fieldSize else { return }
        let range = max(1e-12, result.finiteMaximum - result.finiteMinimum)
        let cellWidth = size.width / CGFloat(result.width)
        let cellHeight = size.height / CGFloat(result.height)
        for row in 0..<result.height {
            for column in 0..<result.width {
                let value = result.values[start + row * result.width + column]
                guard value.isFinite else { continue }
                let normalized = ((value - result.finiteMinimum) / range).clamped(to: 0...1)
                let color = fieldColor(normalized)
                context.fill(
                    Path(CGRect(x: CGFloat(column) * cellWidth, y: CGFloat(row) * cellHeight, width: cellWidth + 0.6, height: cellHeight + 0.6)),
                    with: .color(color)
                )
            }
        }
    }

    private func drawSpatialPath(_ result: SimulationResult, date: Date, context: inout GraphicsContext, size: CGSize) {
        let count = result.values.count / 3
        guard count > 1 else { return }
        let fraction = reduceMotion ? 1 : 0.42 + 0.58 * (0.5 + 0.5 * sin(date.timeIntervalSinceReferenceDate * 0.55))
        let visible = max(2, Int(Double(count) * fraction))
        let strideBy = max(1, visible / 2_800)
        var points = [(Double, Double)]()
        points.reserveCapacity(visible / strideBy)
        for index in Swift.stride(from: 0, to: visible, by: strideBy) {
            let x = result.values[index * 3]
            let z = result.values[index * 3 + 2]
            if x.isFinite, z.isFinite { points.append((x, z)) }
        }
        drawNormalizedPath(points, context: &context, size: size, glow: true)
    }

    private func drawTriangles(_ result: SimulationResult, date: Date, context: inout GraphicsContext, size: CGSize) {
        guard result.values.count >= 18, result.values.first != 0 else { return }
        let triangleCount = (result.values.count - 1) / 18
        let visible = reduceMotion ? triangleCount : max(1, Int(Double(triangleCount) * (0.68 + 0.32 * abs(sin(date.timeIntervalSinceReferenceDate * 0.45)))))
        let strideBy = max(1, triangleCount / 1_500)
        var triangles = [[(Double, Double)]]()
        for triangle in Swift.stride(from: 0, to: min(visible, triangleCount), by: strideBy) {
            let base = 1 + triangle * 18
            var points = [(Double, Double)]()
            for vertex in 0..<3 {
                let offset = base + vertex * 6
                let x = result.values[offset]
                let z = result.values[offset + 2]
                if x.isFinite, z.isFinite { points.append((x, z)) }
            }
            if points.count == 3 { triangles.append(points) }
        }
        let flattened = triangles.flatMap { $0 }
        guard !flattened.isEmpty else { return }
        let minX = flattened.map(\.0).min() ?? 0
        let maxX = flattened.map(\.0).max() ?? 1
        let minY = flattened.map(\.1).min() ?? 0
        let maxY = flattened.map(\.1).max() ?? 1
        let dx = max(1e-9, maxX - minX)
        let dy = max(1e-9, maxY - minY)
        func projected(_ point: (Double, Double)) -> CGPoint {
            CGPoint(
                x: size.width * (0.08 + 0.84 * CGFloat((point.0 - minX) / dx)),
                y: size.height * (0.92 - 0.84 * CGFloat((point.1 - minY) / dy))
            )
        }
        var path = Path()
        for triangle in triangles {
            path.move(to: projected(triangle[0]))
            path.addLine(to: projected(triangle[1]))
            path.addLine(to: projected(triangle[2]))
            path.closeSubpath()
        }
        context.stroke(
            path,
            with: .linearGradient(
                Gradient(colors: [ForgeTheme.cyan, ForgeTheme.violet, ForgeTheme.emerald]),
                startPoint: .zero,
                endPoint: CGPoint(x: size.width, y: size.height)
            ),
            style: StrokeStyle(lineWidth: 0.9, lineJoin: .round)
        )
    }

    private func drawSignal(_ result: SimulationResult, context: inout GraphicsContext, size: CGSize) {
        let finite = result.values.enumerated().filter { $0.element.isFinite }
        guard finite.count > 1 else { return }
        let strideBy = max(1, finite.count / 1_400)
        let range = max(1e-12, result.finiteMaximum - result.finiteMinimum)
        var path = Path()
        var started = false
        for item in finite.enumerated() where item.offset.isMultiple(of: strideBy) {
            let x = CGFloat(item.offset) / CGFloat(max(1, finite.count - 1)) * size.width
            let n = (item.element.element - result.finiteMinimum) / range
            let y = size.height * (0.88 - CGFloat(n.clamped(to: 0...1)) * 0.76)
            if started { path.addLine(to: CGPoint(x: x, y: y)) } else { path.move(to: CGPoint(x: x, y: y)); started = true }
        }
        context.stroke(path, with: .color(accent.opacity(0.18)), style: StrokeStyle(lineWidth: 8, lineCap: .round, lineJoin: .round))
        context.stroke(path, with: .linearGradient(Gradient(colors: [accent, ForgeTheme.violet, ForgeTheme.emerald]), startPoint: .zero, endPoint: CGPoint(x: size.width, y: size.height)), style: StrokeStyle(lineWidth: 2.2, lineCap: .round, lineJoin: .round))
    }

    private func drawCampaign(_ result: SimulationResult, context: inout GraphicsContext, size: CGSize) {
        let finite = result.values.filter(\.isFinite)
        guard !finite.isEmpty else { return }
        let count = min(28, finite.count)
        let range = max(1e-12, result.finiteMaximum - result.finiteMinimum)
        let gap: CGFloat = 5
        let barWidth = max(3, (size.width - CGFloat(count + 1) * gap) / CGFloat(count))
        for index in 0..<count {
            let value = finite[index]
            let n = ((value - result.finiteMinimum) / range).clamped(to: 0...1)
            let height = max(5, CGFloat(n) * size.height * 0.62)
            let rect = CGRect(x: gap + CGFloat(index) * (barWidth + gap), y: size.height * 0.82 - height, width: barWidth, height: height)
            context.fill(Path(roundedRect: rect, cornerRadius: min(5, barWidth / 2)), with: .linearGradient(Gradient(colors: [accent.opacity(0.35), accent]), startPoint: CGPoint(x: rect.midX, y: rect.maxY), endPoint: CGPoint(x: rect.midX, y: rect.minY)))
        }
        let center = CGPoint(x: size.width * 0.75, y: size.height * 0.24)
        for ring in 0..<3 {
            let radius = CGFloat(24 + ring * 20)
            context.stroke(Path(ellipseIn: CGRect(x: center.x - radius, y: center.y - radius, width: radius * 2, height: radius * 2)), with: .color((ring.isMultiple(of: 2) ? ForgeTheme.emerald : ForgeTheme.violet).opacity(0.35)), lineWidth: 1)
        }
    }

    private func drawIdle(in context: inout GraphicsContext, size: CGSize) {
        let center = CGPoint(x: size.width / 2, y: size.height / 2)
        for ring in 0..<5 {
            let radius = CGFloat(26 + ring * 22)
            context.stroke(Path(ellipseIn: CGRect(x: center.x - radius, y: center.y - radius, width: radius * 2, height: radius * 2)), with: .color(accent.opacity(0.10 + Double(ring) * 0.025)), lineWidth: 1)
        }
        context.draw(Text("CHOOSE A KERNEL").font(.system(size: 13, weight: .bold, design: .monospaced)).foregroundColor(accent), at: center)
    }

    private func drawNormalizedPath(_ points: [(Double, Double)], context: inout GraphicsContext, size: CGSize, glow: Bool) {
        guard points.count > 1 else { return }
        let xs = points.map(\.0), ys = points.map(\.1)
        let minX = xs.min() ?? 0, maxX = xs.max() ?? 1, minY = ys.min() ?? 0, maxY = ys.max() ?? 1
        let dx = max(1e-9, maxX - minX), dy = max(1e-9, maxY - minY)
        var path = Path()
        for (index, point) in points.enumerated() {
            let p = CGPoint(x: size.width * (0.08 + 0.84 * CGFloat((point.0 - minX) / dx)), y: size.height * (0.92 - 0.84 * CGFloat((point.1 - minY) / dy)))
            if index == 0 { path.move(to: p) } else { path.addLine(to: p) }
        }
        if glow { context.stroke(path, with: .color(accent.opacity(0.16)), lineWidth: 9) }
        context.stroke(path, with: .linearGradient(Gradient(colors: [ForgeTheme.cyan, ForgeTheme.violet, ForgeTheme.emerald]), startPoint: .zero, endPoint: CGPoint(x: size.width, y: size.height)), style: StrokeStyle(lineWidth: glow ? 2.1 : 1, lineCap: .round, lineJoin: .round))
    }

    private func fieldColor(_ value: Double) -> Color {
        if value < 0.33 {
            return Color(red: 0.02, green: 0.14 + value * 0.8, blue: 0.24 + value * 1.7)
        } else if value < 0.68 {
            let t = (value - 0.33) / 0.35
            return Color(red: 0.08 + 0.34 * t, green: 0.55 + 0.35 * t, blue: 0.78 - 0.32 * t)
        } else {
            let t = (value - 0.68) / 0.32
            return Color(red: 0.42 + 0.55 * t, green: 0.90 - 0.18 * t, blue: 0.46 - 0.22 * t)
        }
    }
}

private struct KernelActivityOverlay: View {
    let accent: Color
    @Environment(\.accessibilityReduceMotion) private var reduceMotion

    var body: some View {
        TimelineView(.animation(minimumInterval: reduceMotion ? 1 : 1 / 30)) { timeline in
            Canvas { context, size in
                let time = reduceMotion ? 0 : timeline.date.timeIntervalSinceReferenceDate
                let center = CGPoint(x: size.width / 2, y: size.height / 2)
                for index in 0..<24 {
                    let phase = time * 1.3 + Double(index) * 0.61
                    let radius = CGFloat(42 + (index % 6) * 16)
                    let point = CGPoint(
                        x: center.x + CGFloat(cos(phase)) * radius,
                        y: center.y + CGFloat(sin(phase * 1.17)) * radius * 0.58
                    )
                    context.fill(Path(ellipseIn: CGRect(x: point.x - 2.2, y: point.y - 2.2, width: 4.4, height: 4.4)), with: .color((index.isMultiple(of: 3) ? ForgeTheme.violet : accent).opacity(0.74)))
                }
            }
        }
        .background(.black.opacity(0.18))
        .overlay(alignment: .bottom) {
            VStack(spacing: 4) {
                Text("SOLVING THE DECLARED KERNEL")
                    .font(.system(size: ForgeTheme.size(12), weight: .bold, design: .monospaced))
                    .foregroundStyle(accent)
                Text("assembling · iterating · checking the result packet")
                    .font(.system(size: ForgeTheme.size(11), design: .rounded))
                    .foregroundStyle(ForgeTheme.secondary)
            }
            .padding(.bottom, 18)
        }
    }
}

private extension Double {
    func clamped(to range: ClosedRange<Double>) -> Double { min(range.upperBound, max(range.lowerBound, self)) }
}
