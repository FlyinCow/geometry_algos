#pragma once
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>

namespace ga {
template <size_t Dim = 3, std::floating_point T = double>
    requires(Dim >= 2)
class Vec {
public:
    using elem_type = T;
    template <std::convertible_to<T>... Args>
        requires(sizeof...(Args) == Dim)
    constexpr Vec(Args... args) noexcept : coord{static_cast<T>(args)...} {
    }

    explicit Vec(std::array<T, Dim> coord) noexcept : coord{coord} {
    }

    constexpr Vec() noexcept : coord{} {
    }

    constexpr Vec &operator=(Vec &&) = default;
    constexpr Vec &operator=(const Vec &) = default;
    constexpr Vec(const Vec &) = default;
    constexpr Vec(Vec &&) = default;

    // accessors

    constexpr auto &&x(this auto &&self) noexcept
        requires(Dim >= 2 && Dim <= 4)
    {
        return self.coord[0];
    }

    constexpr auto &&y(this auto &&self) noexcept
        requires(Dim >= 2 && Dim <= 4)
    {
        return self.coord[1];
    }

    constexpr auto &&u(this auto &&self) noexcept
        requires(Dim == 2)
    {
        return self.coord[0];
    }

    constexpr auto &&v(this auto &&self) noexcept
        requires(Dim == 2)
    {
        return self.coord[1];
    }

    constexpr auto &&z(this auto &&self) noexcept
        requires(Dim == 3 || Dim == 4)
    {
        return self.coord[2];
    }

    constexpr auto &&w(this auto &&self) noexcept
        requires(Dim == 4)
    {
        return self.coord[3];
    }

    // [WARN]: unchecked accsess
    constexpr auto &&operator[](this auto &&self, size_t i) noexcept {
        return self.coord[i];
    }

    // operators

    constexpr Vec &operator+=(const Vec &other) noexcept {
        for (size_t i = 0; i < Dim; i++) {
            coord[i] += other.coord[i];
        }
        return *this;
    }

    constexpr friend Vec operator+(Vec lhs, const Vec &rhs) noexcept {
        lhs += rhs;
        return lhs;
    }

    constexpr Vec &operator*=(T k) noexcept {
        for (size_t i = 0; i < Dim; i++) {
            coord[i] *= k;
        }
        return *this;
    }
    constexpr friend Vec operator*(Vec v, T k) noexcept {
        v *= k;
        return v;
    }
    constexpr friend Vec operator*(T k, Vec v) noexcept {
        v *= k;
        return v;
    }

    constexpr Vec &operator/=(T k) noexcept {
        for (size_t i = 0; i < Dim; i++) {
            coord[i] /= k;
        }
        return *this;
    }
    constexpr friend Vec operator/(Vec v, T k) noexcept {
        v /= k;
        return v;
    }

    constexpr Vec operator-() const noexcept {
        Vec result{};
        for (size_t i = 0; i < Dim; i++) {
            result[i] = -coord[i];
        }
        return result;
    }

    constexpr Vec &operator-=(const Vec &other) noexcept {
        for (size_t i = 0; i < Dim; i++) {
            coord[i] -= other.coord[i];
        }
        return *this;
    }

    constexpr friend Vec operator-(Vec lhs, const Vec &rhs) noexcept {
        lhs -= rhs;
        return lhs;
    }

    constexpr T dot(Vec other) const noexcept {
        T res{};
        for (size_t i = 0; i < Dim; i++) {
            res += coord[i] * other[i];
        }
        return res;
    }

    constexpr T square_norm() const noexcept {
        return dot(*this);
    }

    T norm() const {
        return std::sqrt(square_norm()); // not constexpr noexcept
    }

    constexpr T square_distance_to(this const Vec &self, const Vec &other) noexcept {
        return (self - other).square_norm();
    }

    constexpr T distance_to(this const Vec &self, const Vec &other) noexcept {
        return std::sqrt(self.square_distance_to(other));
    }

    // [WARN]: could overflow for big value
    constexpr Vec normalized(this const Vec &self) {
        return self / self.norm();
    }

    // [WARN]: could overflow for big value
    constexpr Vec &normalize(this Vec &self) {
        self /= self.norm();
        return self;
    }

    constexpr Vec cross(Vec other) const noexcept
        requires(Dim == 3)
    {
        return Vec{
            coord[1] * other[2] - coord[2] * other[1],
            coord[2] * other[0] - coord[0] * other[2],
            coord[0] * other[1] - coord[1] * other[0]};
    }

    constexpr bool is_orthogonal_to(this const Vec &self, const Vec &other, T eps) noexcept {
        return std::abs(self.dot(other)) <= eps;
    }

    constexpr bool is_zero(T eps) {
        return this->square_distance_to({}) <= eps * eps;
    }

private:
    std::array<T, Dim> coord;
};

using Vec3d = Vec<>;
using Vec3f = Vec<3, float>;
using Vec2d = Vec<2>;
using Vec2f = Vec<2, float>;

template <size_t Dim, class T>
T squared_distance(const Vec<Dim, T> &a, const Vec<Dim, T> &b) noexcept {
    return (b - a).square_norm();
}

template <size_t Dim, class T>
T distance(const Vec<Dim, T> &a, const Vec<Dim, T> &b) noexcept {
    return (b - a).norm();
}

template <size_t Dim, class T>
bool is_orthogonal(const Vec<Dim, T> &a, const Vec<Dim, T> &b, T eps) noexcept {
    return std::abs(a.dot(b)) <= eps;
    // return a.is_orthogonal_to(b, eps);
}

// [TODO]: linear relative

} // namespace ga

template <size_t Dim, class T>
struct std::formatter<ga::Vec<Dim, T>> {
    enum class BracketType { None,
                             Parenthesis,
                             Square,
                             Curly };
    enum class SepType { Space,
                         Comma };

    BracketType bracket = BracketType::Parenthesis;
    SepType sep = SepType::Space;
    int precision = 8;

    constexpr auto parse(std::format_parse_context &ctx) {
        auto pos = ctx.begin();
        while (pos != ctx.end() && *pos != '}') {
            switch (*pos) {
            case 'p':
                bracket = BracketType::Parenthesis;
                break;
            case 's':
                bracket = BracketType::Square;
                break;
            case 'c':
                bracket = BracketType::Curly;
                break;
            case 'n':
                bracket = BracketType::None;
                break;
            case ',':
                sep = SepType::Comma;
                break;
            case 'g': {
                ++pos;
                if (pos != ctx.end() && std::isdigit(*pos)) {
                    int p = *pos - '0';
                    if (p > 0)
                        precision = p;
                }
                break;
            }
            }
            ++pos;
        }
        return pos;
    }

    auto format(const ga::Vec<Dim, T> &v, std::format_context &ctx) const {
        std::format_context::iterator out = ctx.out();

        switch (bracket) {
        case BracketType::Parenthesis: out = std::format_to(out, "(");  break;
        case BracketType::Square:      out = std::format_to(out, "[");  break;
        case BracketType::Curly:       out = std::format_to(out, "{{"); break;
        case BracketType::None:        break;
        }

        const std::string_view separator = sep == SepType::Comma ? ", " : " ";
        for (size_t i = 0; i < Dim; ++i) {
            out = std::format_to(out, "{:.{}g}", v[i], precision);
            if (i + 1 < Dim)
                out = std::format_to(out, "{}", separator);
        }

        switch (bracket) {
        case BracketType::Parenthesis: out = std::format_to(out, ")");  break;
        case BracketType::Square:      out = std::format_to(out, "]");  break;
        case BracketType::Curly:       out = std::format_to(out, "}}"); break;
        case BracketType::None:        break;
        }

        return out;
    }
};