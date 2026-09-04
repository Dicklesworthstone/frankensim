import SwiftUI
import UIKit

enum ForgeAppearance: String {
    static let storageKey = "frankensim.appearance"
    case dark
    case light
    var colorScheme: ColorScheme { self == .dark ? .dark : .light }
}

enum ForgeTheme {
    static let textScaleStorageKey = "frankensim.textScale"
    static let defaultTextScale = 1.0
    static let minimumTextScale = 0.8
    static let maximumTextScale = 1.5
    static let textScaleStep = 0.1

    static let background = adaptive(dark: UIColor(red: 0.012, green: 0.024, blue: 0.032, alpha: 1), light: UIColor(red: 0.93, green: 0.96, blue: 0.972, alpha: 1))
    static let panel = adaptive(dark: UIColor(red: 0.025, green: 0.050, blue: 0.062, alpha: 1), light: UIColor(red: 0.985, green: 0.995, blue: 1, alpha: 1))
    static let raised = adaptive(dark: UIColor(red: 0.038, green: 0.073, blue: 0.086, alpha: 1), light: UIColor(red: 0.84, green: 0.91, blue: 0.935, alpha: 1))
    static let stroke = adaptive(dark: UIColor(white: 1, alpha: 0.06), light: UIColor(red: 0.03, green: 0.22, blue: 0.29, alpha: 0.16))
    static let text = adaptive(dark: UIColor(red: 0.93, green: 0.96, blue: 0.98, alpha: 1), light: UIColor(red: 0.035, green: 0.09, blue: 0.125, alpha: 1))
    static let secondary = adaptive(dark: UIColor(red: 0.62, green: 0.70, blue: 0.76, alpha: 1), light: UIColor(red: 0.27, green: 0.35, blue: 0.39, alpha: 1))
    static let emerald = adaptive(dark: UIColor(red: 0.25, green: 0.90, blue: 0.65, alpha: 1), light: UIColor(red: 0.015, green: 0.415, blue: 0.255, alpha: 1))
    static let cyan = adaptive(dark: UIColor(red: 0.22, green: 0.80, blue: 0.96, alpha: 1), light: UIColor(red: 0.015, green: 0.39, blue: 0.545, alpha: 1))
    static let violet = adaptive(dark: UIColor(red: 0.64, green: 0.46, blue: 0.98, alpha: 1), light: UIColor(red: 0.39, green: 0.22, blue: 0.68, alpha: 1))
    static let amber = adaptive(dark: UIColor(red: 0.98, green: 0.69, blue: 0.25, alpha: 1), light: UIColor(red: 0.66, green: 0.35, blue: 0.01, alpha: 1))
    static let coral = adaptive(dark: UIColor(red: 0.97, green: 0.42, blue: 0.48, alpha: 1), light: UIColor(red: 0.70, green: 0.12, blue: 0.18, alpha: 1))

    private static func adaptive(dark: UIColor, light: UIColor) -> Color {
        Color(uiColor: UIColor { traits in traits.userInterfaceStyle == .dark ? dark : light })
    }

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
        let textScale = CGFloat(storedTextScale)
#if targetEnvironment(macCatalyst)
        return base * 1.22 * textScale
#else
        return UIFontMetrics(forTextStyle: .body).scaledValue(for: base) * textScale
#endif
    }

    static var storedTextScale: Double {
        let defaults = UserDefaults.standard
        guard defaults.object(forKey: textScaleStorageKey) != nil else { return defaultTextScale }
        return normalizedTextScale(defaults.double(forKey: textScaleStorageKey))
    }

    static func normalizedTextScale(_ candidate: Double) -> Double {
        guard candidate.isFinite else { return defaultTextScale }
        let clamped = min(max(candidate, minimumTextScale), maximumTextScale)
        return (clamped / textScaleStep).rounded() * textScaleStep
    }

    static func adjustedTextScale(from current: Double, steps: Int) -> Double {
        normalizedTextScale(current + Double(steps) * textScaleStep)
    }

    static func dynamicTypeSize(for scale: Double) -> DynamicTypeSize {
        switch normalizedTextScale(scale) {
        case ..<0.9: return .small
        case ..<1.0: return .medium
        case ..<1.1: return .large
        case ..<1.2: return .xLarge
        case ..<1.3: return .xxLarge
        case ..<1.4: return .xxxLarge
        default: return .accessibility1
        }
    }
}

struct ForgeAppearanceButton: View {
    @Binding var selection: String
    private var appearance: ForgeAppearance { ForgeAppearance(rawValue: selection) ?? .dark }

    var body: some View {
        Button {
            selection = appearance == .dark ? ForgeAppearance.light.rawValue : ForgeAppearance.dark.rawValue
        } label: {
            Image(systemName: appearance == .dark ? "sun.max.fill" : "moon.stars.fill")
                .frame(width: 44, height: 44)
        }
        .foregroundStyle(appearance == .dark ? ForgeTheme.amber : ForgeTheme.cyan)
        .accessibilityIdentifier("appearance-toggle")
        .accessibilityLabel(appearance == .dark ? "Switch to light mode" : "Switch to dark mode")
        .accessibilityValue(appearance == .dark ? "Dark mode" : "Light mode")
        .accessibilityHint("Remembers this choice for future launches")
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
                            colors: [accent.opacity(0.42), ForgeTheme.stroke, ForgeTheme.violet.opacity(0.18)],
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
