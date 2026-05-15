#include "math/Vec.h"
#include <print>
#include <type_traits>

using namespace ga;

struct Opened {};
struct Connected {};
struct Closed {};

template <class State>
struct Connection {
    Connection<Opened> open()
        requires(std::is_same_v<State, Closed>);
    Connection<Connected> connect()
        requires(std::is_same_v<State, Opened>);
};

int main() {
    auto v = Vec2d{3, 4};
}