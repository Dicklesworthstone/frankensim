import SwiftUI
import UIKit

enum ForgeTheme {
    static let background = Color(red: 0.012, green: 0.024, blue: 0.032)
    static let panel = Color(red: 0.025, green: 0.050, blue: 0.062)
    static let raised = Color(red: 0.038, green: 0.073, blue: 0.086)
    static let text = Color(red: 0.93, green: 0.96, blue: 0.98)
    static let secondary = Color(red: 0.62, green: 0.70, blue: 0.76)
    static let emerald = Color(red: 0.25, green: 0.90, blue: 0.65)
    static let cyan = Color(red: 0.22, green: 0.80, blue: 0.96)
    static let violet = Color(red: 0.64, green: 0.46, blue: 0.98)
    static let amber = Color(red: 0.98, green: 0.69, blue: 0.25)
    static let coral = Color(red: 0.97, green: 0.42, blue: 0.48)

    static func accent(_ family: AccentFamily) -> Color {
        switch family {
        case .cyan: cyan
        case .violet: violet
        case .emerald: emerald
        case .amber: amber
        case .coral: coral
        }
    }

    static func size(_ base: CGFloat) -> CGFloat {
#if targetEnvironment(macCatalyst)
        base * 1.22
#else
        UIFontMetrics(forTextStyle: .body).scaledValue(for: base)
#endif
    }
}

enum AccentFamily: String, Sendable {
    case cyan, violet, emerald, amber, coral
}

struct ForgeBackground: View {
    @Environment(\.accessibilityReduceTransparency) private var reduceTransparency

    var body: some View {
        ZStack {
            ForgeTheme.background
            RadialGradient(
                colors: [ForgeTheme.emerald.opacity(reduceTransparency ? 0.05 : 0.15), .clear],
                center: .topLeading,
                startRadius: 0,
                endRadius: 760
            )
            RadialGradient(
                colors: [ForgeTheme.violet.opacity(reduceTransparency ? 0.03 : 0.10), .clear],
                center: .bottomTrailing,
                startRadius: 0,
                endRadius: 850
            )
            Canvas { context, size in
                let step: CGFloat = 44
                var path = Path()
                stride(from: CGFloat.zero, through: size.width, by: step).forEach { x in
                    path.move(to: CGPoint(x: x, y: 0))
                    path.addLine(to: CGPoint(x: x, y: size.height))
                }
                stride(from: CGFloat.zero, through: size.height, by: step).forEach { y in
                    path.move(to: CGPoint(x: 0, y: y))
                    path.addLine(to: CGPoint(x: size.width, y: y))
                }
                context.stroke(path, with: .color(ForgeTheme.cyan.opacity(0.028)), lineWidth: 0.6)
            }
        }
        .ignoresSafeArea()
        .accessibilityHidden(true)
    }
}

struct ForgePanel<Content: View>: View {
    let accent: Color
    let content: Content

    init(accent: Color = ForgeTheme.cyan, @ViewBuilder content: () -> Content) {
        self.accent = accent
        self.content = content()
    }

    var body: some View {
        content
            .padding(16)
            .background(ForgeTheme.panel.opacity(0.92), in: RoundedRectangle(cornerRadius: 22, style: .continuous))
            .overlay {
                RoundedRectangle(cornerRadius: 22, style: .continuous)
                    .stroke(
                        LinearGradient(
                            colors: [accent.opacity(0.42), Color.white.opacity(0.06), ForgeTheme.violet.opacity(0.18)],
                            startPoint: .topLeading,
                            endPoint: .bottomTrailing
                        ),
                        lineWidth: 1
                    )
            }
            .shadow(color: accent.opacity(0.08), radius: 22, y: 10)
    }
}

struct FrankenSimWordmark: View {
    var body: some View {
        (
            Text("F")
                .font(.system(size: ForgeTheme.size(23), weight: .black, design: .rounded))
                .foregroundColor(ForgeTheme.text)
            + Text("RANKEN")
                .font(.system(size: ForgeTheme.size(15), weight: .black, design: .rounded))
                .foregroundColor(ForgeTheme.text)
            + Text("S")
                .font(.system(size: ForgeTheme.size(23), weight: .black, design: .rounded))
                .foregroundColor(ForgeTheme.emerald)
            + Text("IM")
                .font(.system(size: ForgeTheme.size(15), weight: .black, design: .rounded))
                .foregroundColor(ForgeTheme.emerald)
        )
        .kerning(0.7)
        .lineLimit(1)
        .minimumScaleFactor(0.68)
        .accessibilityElement(children: .ignore)
        .accessibilityLabel("FrankenSim")
    }
}

struct PrimaryForgeButtonStyle: ButtonStyle {
    let tint: Color
    @Environment(\.isEnabled) private var isEnabled

    func makeBody(configuration: Configuration) -> some View {
        configuration.label
            .font(.system(size: ForgeTheme.size(14), weight: .bold, design: .rounded))
            .foregroundStyle(ForgeTheme.background)
            .padding(.horizontal, 18)
            .frame(minHeight: 48)
            .background(tint.opacity(configuration.isPressed ? 0.78 : 1), in: Capsule())
            .scaleEffect(configuration.isPressed ? 0.97 : 1)
            .opacity(isEnabled ? 1 : 0.38)
            .animation(.easeOut(duration: 0.14), value: configuration.isPressed)
    }
}

struct SecondaryForgeButtonStyle: ButtonStyle {
    let tint: Color

    func makeBody(configuration: Configuration) -> some View {
        configuration.label
            .font(.system(size: ForgeTheme.size(13), weight: .semibold, design: .rounded))
            .foregroundStyle(tint)
            .padding(.horizontal, 14)
            .frame(minHeight: 44)
            .background(tint.opacity(configuration.isPressed ? 0.17 : 0.07), in: Capsule())
            .overlay(Capsule().stroke(tint.opacity(0.35), lineWidth: 1))
    }
}

struct CatalystWindowFreedom: UIViewControllerRepresentable {
    func makeUIViewController(context: Context) -> Controller { Controller() }
    func updateUIViewController(_ controller: Controller, context: Context) { controller.configure() }

    final class Controller: UIViewController {
        override func viewDidAppear(_ animated: Bool) {
            super.viewDidAppear(animated)
            configure()
        }

        override func viewDidLayoutSubviews() {
            super.viewDidLayoutSubviews()
            configure()
        }

        func configure() {
#if targetEnvironment(macCatalyst)
            guard let restrictions = view.window?.windowScene?.sizeRestrictions else { return }
            restrictions.minimumSize = CGSize(width: 760, height: 560)
            restrictions.maximumSize = CGSize(width: 10_000, height: 10_000)
#endif
        }
    }
}
